#!/usr/bin/env python3
"""
The peak check (skyline_peaks.py, SKILL.md step 8f; references/skyline.md).

A feasibility trial on HIVE (2026-10-20) ran Skyline's command line on two finished DIA-NN searches
and found the rules these tests hold the code to:
  * Skyline's tolerance must be DIA-NN's own mass accuracy for the search -- a guessed 10 ppm gave
    gappy, broken peaks on a run DIA-NN had measured at 24 ppm;
  * the mass analyser is always set WITH its resolution (DIA alone left the product analyser at
    Skyline's default "qit" 0.7);
  * one raw file per SLURM array task (--import-no-join): two readers on one timsTOF .d failed with
    "database is locked"; and never --import-process-count, which hung under Wine;
  * Skyline reuses DIA-NN's peak edges by default, so a matching peak is not independent
    confirmation: the report says so and never "confirmed by Skyline"; the warning signs are the
    precursor-runs where Skyline picked a different peak (apex > 6 s away, or outside DIA-NN's
    edges, or no peak);
  * a step that cannot run (no container) says so in the record, the report and MANIFEST.txt.

Stdlib only; the pictures test needs matplotlib and is skipped without it.
"""
import csv
import io
import json
import ntpath
import os
import stat
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import client_text      # noqa: E402
import scratch_files    # noqa: E402
import skyline_peaks as sp  # noqa: E402
from resolve_defaults import SKYLINE_PIN  # noqa: E402

try:
    import matplotlib  # noqa: F401  -- only whether it is there
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False

CMD = ("/opt/diann/diann-linux --qvalue 0.01 --xic 10 --cut K*,R* --missed-cleavages 1 "
       "--min-pep-len 7 --max-pep-len 30 --min-pr-mz 357 --max-pr-mz 1105 --min-pr-charge 2 "
       "--max-pr-charge 4 --min-fr-mz 200 --max-fr-mz 1800 --met-excision --unimod4 "
       "--f /data/run A/r1.raw --f /data/r2.raw --fasta /db/search.fasta --out /s/report.parquet")
AUTO_LOG = ("DIA-NN 2.7.0 Academia  (Data-Independent Acquisition by Neural Networks)\n"
            "Compiled on Sep 15 2026\nLogical CPU cores: 16\n" + CMD + "\n\n"
            "DIA-NN will automatically optimise the mass accuracy for the first run of the "
            "experiment\n[0:57] Recommended MS1 mass accuracy setting: 8 ppm\n"
            "[1:45] Optimised mass accuracy: 24 ppm\n[4:20] Recommended MS1 mass accuracy "
            "setting: 7 ppm\n")
FIXED_LOG = ("DIA-NN 2.7.0 Academia\n" + CMD.replace("--unimod4", "--unimod4 --mass-acc 15 "
                                                    "--mass-acc-ms1 12") + "\n"
             "Mass accuracy will be fixed to 1.5e-05 (MS2) and 1.2e-05 (MS1)\n"
             "[1:25] Recommended MS1 mass accuracy setting: 12 ppm\n")


def plan_for(family="orbitrap", raw=None, ms1=24.0, ms2=24.0, work="/w/skyline_work"):
    raw = raw or ["/data/run A/r1.raw", "/data/r2.raw", "/data/r3.raw"]
    out = os.path.dirname(work)
    return {"work": work, "out": out, "doc": os.path.join(work, "T_top_hits.sky"),
            "share_zip": os.path.join(out, "T_top_hits.sky.zip"), "raw_files": raw,
            "family": family, "mass_accuracy": {"ms1_ppm": ms1, "ms2_ppm": ms2, "source": "test"},
            "settings": sp.search_settings(sp.command_flags(CMD)),
            "speclib_link": os.path.join(work, "search", "report-lib.parquet.skyline.speclib"),
            "targets_fasta": os.path.join(work, "targets.fasta"), "sif": "/c/sky.sif",
            "ro_binds": ["/data", "/s"], "python": "/env/bin/python",
            "script": "/skill/scripts/skyline_peaks.py", "max_parallel": 8, "name": "T"}


def opts(args):
    """The value of each --flag=value in a SkylineCmd argument list."""
    return {k: v for k, _, v in (a.partition("=") for a in args)}


class TestMassAccuracy(unittest.TestCase):
    """DIA-NN's mass accuracy for THIS search, never a fixed number."""

    def test_measured_by_diann(self):
        ma = sp.mass_accuracy(AUTO_LOG)
        self.assertEqual((ma["ms1_ppm"], ma["ms2_ppm"]), (24.0, 24.0))
        self.assertIn("Optimised mass accuracy: 24 ppm", ma["source"])
        self.assertIn("recommendation for MS1 (7, 8 ppm)", ma["source"])

    def test_fixed_in_the_search(self):
        ma = sp.mass_accuracy(FIXED_LOG)
        self.assertEqual((ma["ms1_ppm"], ma["ms2_ppm"]), (12.0, 15.0))
        self.assertIn("fixed in the DIA-NN search", ma["source"])

    def test_pinned_in_the_command_without_the_fixed_line(self):
        log = "DIA-NN 2.7.0\n" + CMD + " --mass-acc 18\n"
        ma = sp.mass_accuracy(log)
        self.assertEqual((ma["ms1_ppm"], ma["ms2_ppm"]), (18.0, 18.0))
        self.assertIn("--mass-acc", ma["source"])

    def test_search_provenance_as_the_last_record(self):
        prov = {"mass_acc": {"mode": "measured", "ms1": 9.0, "ms2": 14.0, "reason": "probe"}}
        ma = sp.mass_accuracy("DIA-NN 2.7.0\n" + CMD + "\n", prov)
        self.assertEqual((ma["ms1_ppm"], ma["ms2_ppm"]), (9.0, 14.0))
        self.assertIn("search_provenance.json", ma["source"])

    def test_nothing_recorded_is_not_run_never_a_default(self):
        with self.assertRaises(sp.Skip) as e:
            sp.mass_accuracy("DIA-NN 2.7.0\n" + CMD + "\n", {})
        self.assertIn("never" if "never" in str(e.exception) else "10 ppm", str(e.exception))
        self.assertIn("--ppm-ms2", str(e.exception))

    def test_the_users_value_needs_both_levels(self):
        self.assertEqual(sp.mass_accuracy(AUTO_LOG, None, 5, 9)["ms2_ppm"], 9.0)
        with self.assertRaises(sp.Skip):
            sp.mass_accuracy(AUTO_LOG, None, None, 9)


class TestSearchSettings(unittest.TestCase):
    def test_copied_from_the_diann_command(self):
        s = sp.search_settings(sp.command_flags(CMD))
        self.assertEqual(s["enzyme"], "Trypsin/P")
        self.assertEqual((s["missed_cleavages"], s["min_length"], s["max_length"]), (1, 7, 30))
        self.assertEqual(s["charges"], [2, 3, 4])
        self.assertEqual(s["fragment_mz"], [200.0, 1800.0])
        self.assertTrue(any("methionine" in x for x in s["not_carried_over"]))

    def test_a_path_with_a_space_stays_whole(self):
        self.assertIn("/data/run A/r1.raw", sp.command_flags(CMD)["f"])

    def test_enzymes(self):
        self.assertEqual(sp.enzyme_for("K*,R*"), "Trypsin/P")
        self.assertEqual(sp.enzyme_for("K*,R*,!*P"), "Trypsin")
        self.assertEqual(sp.enzyme_for("K*,!*P"), "LysC")
        self.assertEqual(sp.enzyme_for("K*"), "LysC/P")
        self.assertEqual(sp.enzyme_for("*D"), "AspN")
        for bad in ("K*,*D", "X+", "W*,!*P"):
            with self.assertRaises(sp.Skip):
                sp.enzyme_for(bad)

    def test_a_flag_not_stated_is_not_guessed(self):
        with self.assertRaises(sp.Skip) as e:
            sp.search_settings(sp.command_flags(CMD.replace("--missed-cleavages 1 ", "")))
        self.assertIn("--missed-cleavages", str(e.exception))
        with self.assertRaises(sp.Skip):
            sp.search_settings(sp.command_flags(CMD.replace("--unimod4 ", "")))

    def test_instrument_family(self):
        self.assertEqual(sp.instrument_family({"instrument_class": "timstof"}, ["/a/x.d"])[0],
                         "timstof")
        self.assertEqual(sp.instrument_family({"instrument_class": "orbitrap_astral"},
                                              ["/a/x.raw"])[0], "orbitrap")
        self.assertEqual(sp.instrument_family(None, ["/a/x.raw", "/a/y.raw"])[0], "orbitrap")
        for man, raw in (({"instrument_class": "sciex_tof"}, ["/a/x.wiff"]),
                         (None, ["/a/x.mzML"]), ({"instrument_class": "timstof"}, ["/a/x.raw"]),
                         ({"instrument_class": "orbitrap_generic",
                           "resolution": {"ms2_analyzer": "ITMS"}}, ["/a/x.raw"])):
            with self.assertRaises(sp.Skip):
                sp.instrument_family(man, raw)


class TestCommandBuilders(unittest.TestCase):
    def test_orbitrap_centroided_at_diann_ppm(self):
        o = opts(sp.build_args(plan_for(), 2.0))
        self.assertEqual(o["--full-scan-precursor-analyzer"], "centroided")
        self.assertEqual(o["--full-scan-precursor-res"], "24")
        self.assertEqual(o["--full-scan-product-analyzer"], "centroided")
        self.assertEqual(o["--full-scan-product-res"], "24")
        self.assertNotIn("10", (o["--full-scan-precursor-res"], o["--full-scan-product-res"]))

    def test_timstof_tof_resolving_power_from_ppm(self):
        o = opts(sp.build_args(plan_for("timstof", ["/d/a.d", "/d/b.d"], 15, 15), 1.0))
        self.assertEqual(o["--full-scan-precursor-analyzer"], "tof")
        self.assertEqual(o["--full-scan-product-analyzer"], "tof")
        # Skyline's TOF window is +- m/z / R: +-15 ppm is R = 66,667
        self.assertEqual(o["--full-scan-product-res"], "66667")
        self.assertEqual(o["--ims-library-res"], str(sp.IMS_LIBRARY_RES))

    def test_analyser_never_without_its_resolution(self):
        for fam in ("orbitrap", "timstof"):
            o = opts(sp.build_args(plan_for(fam, ["/d/a.d"] if fam == "timstof" else None), 1.0))
            for lvl in ("precursor", "product"):
                self.assertIn(f"--full-scan-{lvl}-analyzer", o)
                self.assertIn(f"--full-scan-{lvl}-res", o)

    def test_one_build_call_with_the_searchs_settings(self):
        p = plan_for()
        args = sp.build_args(p, 2.0)
        o = opts(args)
        self.assertEqual(o["--full-scan-acquisition-method"], "DIA")
        self.assertEqual(o["--pep-digest-enzyme"], "Trypsin/P")
        self.assertEqual(o["--tran-precursor-ion-charges"], "2,3,4")
        self.assertEqual(o["--instrument-min-mz"], "200")
        self.assertEqual(o["--full-scan-rt-filter-tolerance"], "2")
        self.assertTrue(o["--import-search-file"].endswith("report-lib.parquet.skyline.speclib"))
        self.assertTrue(o["--full-scan-isolation-scheme"].startswith("Z:\\data\\run A\\"))
        self.assertEqual(args[-1], "--save")

    def test_windows_paths(self):
        self.assertEqual(sp.win("/quobyte/x y/doc.sky"), "Z:\\quobyte\\x y\\doc.sky")

    def test_one_task_per_raw_file_without_join(self):
        p = plan_for()
        docs = set()
        for i, raw in enumerate(p["raw_files"], 1):
            a = sp.import_args(p, i)
            self.assertIn(f"--import-file={sp.win(raw)}", a)
            self.assertEqual(sum(x.startswith("--import-file") for x in a), 1)
            self.assertIn("--import-no-join", a)
            self.assertNotIn("--save", a)
            docs.add(ntpath.dirname(opts(a)["--in"]))
        # each task its own document folder: Skyline deletes <document folder>/xic/ when an
        # import ends, which killed the other tasks' imports when they shared one
        self.assertEqual(len(docs), len(p["raw_files"]))

    def test_the_import_task_moves_its_cache_beside_the_document(self):
        with tempfile.TemporaryDirectory() as td:
            w = os.path.join(td, "skyline_work")
            os.makedirs(os.path.join(w, "status"))
            os.makedirs(os.path.join(w, "logs"))
            p = dict(plan_for(work=w), out=td)
            with open(p["doc"], "w") as fh:
                fh.write("<srm_settings/>")
            with open(os.path.join(w, "plan.json"), "w") as fh:
                json.dump(p, fh)

            def fake(plan, args, log_name, minutes, rw=()):
                doc = args[0].split("=", 1)[1][2:].replace("\\", "/")     # Z:\a\b -> /a/b
                self.assertTrue(os.path.isfile(doc))                 # the task's copy
                self.assertNotEqual(os.path.dirname(doc), w)
                open(os.path.join(os.path.dirname(doc), "r2_123.raw.skyd"), "w").close()
                return 0, os.path.join(w, "logs", log_name)
            with mock.patch.object(sp, "run_skyline", fake):
                rc = sp.main(["import", "--work", w, "--task", "2"])
            self.assertEqual(rc, 0)
            self.assertTrue(os.path.isfile(os.path.join(w, "r2_123.raw.skyd")))
            with open(os.path.join(w, "status", "import_2.json")) as fh:
                self.assertEqual(json.load(fh)["rc"], 0)
            self.assertFalse(os.path.isdir(os.path.join(w, "tasks", "2")))

    def test_no_import_process_count_anywhere(self):
        p = plan_for()
        queue = ("high", "acct", "acct-high-qos")
        texts = [" ".join(sp.build_args(p, 1.0)), sp.join_batch(p, p["raw_files"]),
                 sp.submit_script(p)] + [" ".join(sp.import_args(p, i)) for i in (1, 2, 3)] \
            + [sp.job_script(p, s, queue) for s in ("build", "import", "join", "prm")]
        for t in texts:
            self.assertNotIn("import-process-count", t)
        with self.assertRaises(AssertionError):
            sp.run_skyline(p, ["--import-process-count=2"], "x.txt", 5)

    def test_join_batch_lists_only_the_imported_files(self):
        p = plan_for()
        batch = sp.join_batch(p, p["raw_files"][:2])
        first = batch.splitlines()[0]
        self.assertEqual(first.count("--import-file="), 2)
        self.assertNotIn("r3.raw", batch)
        self.assertIn('--in="Z:\\w\\skyline_work\\T_top_hits.sky"', first)
        self.assertIn("--share-type=complete", batch)
        self.assertIn('--report-name="' + sp.REPORT_PRECURSORS + '"', batch)
        self.assertIn("--chromatogram-products", batch)
        with self.assertRaises(sp.Skip):
            sp.batch_line(['--in=/a"b'])

    def test_apptainer_binds_raw_and_search_read_only(self):
        argv = sp.apptainer_argv(plan_for(), rw=["/followup"])
        j = " ".join(argv)
        self.assertIn("--bind /w:/w", j)
        self.assertIn("--bind /followup:/followup", j)
        self.assertIn("--bind /data:/data:ro", j)
        self.assertIn("--bind /s:/s:ro", j)
        self.assertEqual(argv[-1], "/c/sky.sif")

        self.assertEqual(sp.ro_binds("/s/search", ["/raw/a/x.d", "/raw/a/y.d", "/raw/b/z.d"],
                                     "/s/search/db.fasta"), ["/raw/a", "/raw/b", "/s/search"])

    def test_job_scripts_follow_ucd_hpc(self):
        p = plan_for()
        queue = ("high", "genome-center-grp", "genome-center-grp-high-qos")
        build = sp.job_script(p, "build", queue)
        arr = sp.job_script(p, "import", queue, array="1-3%3")
        for t in (build, arr):
            for need in ("#SBATCH --account=genome-center-grp", "#SBATCH --partition=high",
                         "#SBATCH --qos=genome-center-grp-high-qos", "#SBATCH --time=",
                         "#SBATCH --mem=16G", "#SBATCH --ntasks=1", "#SBATCH --cpus-per-task="):
                self.assertIn(need, t)
            self.assertNotIn("--exclusive", t)
            self.assertLess(t.index("source /etc/profile.d/modules.sh"), t.index("module load"))
        self.assertIn("_%j.log", build)
        self.assertIn("#SBATCH --array=1-3%3", arr)
        self.assertIn("_%A_%a.log", arr)
        self.assertIn('--task "$SLURM_ARRAY_TASK_ID"', arr)
        low = sp.job_script(p, "join", ("low", "publicgrp", "publicgrp-low-qos"))
        self.assertIn("#SBATCH --requeue", low)

    def test_every_job_runs_a_subcommand_that_exists(self):
        p = plan_for()
        queue = ("high", "acct", "acct-high-qos")
        for step, sub in sp.JOB_SUBCOMMAND.items():
            t = sp.job_script(p, step, queue)
            self.assertIn(f"skyline_peaks.py {sub} --work ", t)
            with mock.patch("sys.stdout", new=io.StringIO()), self.assertRaises(SystemExit) as e:
                sp.main([sub, "--help"])
            self.assertEqual(e.exception.code, 0, sub)
        self.assertIn("prm-run --work", sp.job_script(p, "prm", queue))

    def test_submit_chains_the_three_jobs(self):
        t = sp.submit_script(plan_for())
        self.assertIn("--dependency=afterok:$j1 --kill-on-invalid-dep=yes", t)
        self.assertIn("--dependency=afterany:$j2", t)

    def test_rt_window(self):
        self.assertEqual(sp.rt_window(11.0)[0], 1.0)
        self.assertEqual(sp.rt_window(35.0)[0], 2.0)


def write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def de_row(pg, gene, adj, p, lfc=1.0):
    return {"Protein.Group": pg, "Genes": gene, "logFC": lfc, "P.Value": p, "adj.P.Val": adj,
            "Evidence": "measured in both"}


class TestTopN(unittest.TestCase):
    """ONE rule: the significant proteins by their smallest adjusted p in any comparison, ties by
    the raw p, contaminants out; none significant -> the most abundant measured in every run."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.de = self._td.name
        with open(os.path.join(self.de, "de_provenance.json"), "w") as fh:
            json.dump({"adjp": 0.05, "contrasts": ["A-B", "C-B"]}, fh)

    def tearDown(self):
        self._td.cleanup()

    def plan(self, n=3):
        return {"top": {"value": n, "source": "skill default"}, "de_dir": self.de}

    def test_best_across_contrasts(self):
        write_csv(os.path.join(self.de, "DE_dpc_A.B.csv"), [
            de_row("P1", "G1", 0.01, 0.001), de_row("P2", "G2", 0.001, 1e-5),
            de_row("Cont_P00761", "TRYP", 1e-9, 1e-12), de_row("P3", "G3", 0.2, 0.01),
            de_row("P5", "G5", "", "")])
        write_csv(os.path.join(self.de, "DE_dpc_C.B.csv"), [
            de_row("P1", "G1", 0.0005, 1e-6, -2), de_row("P4", "G4", 0.001, 2e-5),
            de_row("P6;Cont_Q1", "G6", 1e-6, 1e-8)])
        sel = sp.select_targets(self.plan(3))
        self.assertEqual(sel["mode"], "de_hits")
        got = [(t["protein_group"], t["contrast"]) for t in sel["targets"]]
        # P1 by its best comparison (C - B); P2 before P4 on the tie (raw p 1e-5 < 2e-5);
        # contaminant groups and P3 (not significant) left out; N caps the list
        self.assertEqual(got, [("P1", "C vs B"), ("P2", "A vs B"), ("P4", "C vs B")])
        self.assertEqual(sel["adjp"], 0.05)

    def test_no_hits_falls_back_to_abundance_measured_everywhere(self):
        write_csv(os.path.join(self.de, "DE_dpc_A.B.csv"), [de_row("P1", "G1", 0.5, 0.1)])
        write_csv(os.path.join(self.de, "Expression_Matrix.csv"), [
            {"Protein.Group": "P1", "Genes": "G1", "Protein.Names": "", "s1": 20, "s2": 22},
            {"Protein.Group": "P2", "Genes": "G2", "Protein.Names": "", "s1": 30, "s2": 30},
            {"Protein.Group": "P3", "Genes": "G3", "Protein.Names": "", "s1": 25, "s2": 25},
            {"Protein.Group": "Cont_X", "Genes": "K", "Protein.Names": "", "s1": 40, "s2": 40}])
        write_csv(os.path.join(self.de, "Detection_Matrix.csv"), [
            {"Protein.Group": "P1", "s1": 2, "s2": 3}, {"Protein.Group": "P2", "s1": 0, "s2": 4},
            {"Protein.Group": "P3", "s1": 5, "s2": 5}, {"Protein.Group": "Cont_X", "s1": 9, "s2": 9}])
        sel = sp.select_targets(self.plan(5))
        self.assertEqual(sel["mode"], "abundance_no_hits")
        # P2 is the most abundant but inferred in s1: not measured in every run
        self.assertEqual([t["protein_group"] for t in sel["targets"]], ["P3", "P1"])


SKY_HEADER = ["ProteinName", "PeptideSequence", "PeptideModifiedSequence", "PrecursorCharge",
              "IsDecoy", "FileName", "BestRetentionTime", "MinStartTime", "MaxEndTime",
              "TotalArea", "TotalAreaFragment", "LibraryDotProduct", "IsotopeDotProduct"]


def sky_row(prot, mod, z, run, apex, start, end):
    seq = "".join(c for c in mod if c.isalpha() and c.isupper())
    return dict(zip(SKY_HEADER, [prot, seq, mod, z, "False", run + ".raw", apex, start, end,
                                 "1000", "800", "0.9", "0.95"]))


def dia_row(group, mod, z, run, rt, start, stop):
    stripped = "".join(c for c in mod.split("(")[0] if c.isalpha()) if "(" not in mod else \
        "".join(c for i, c in enumerate(mod) if c.isupper() and "(" not in mod[max(0, i - 7):i])
    return {"Run": run, "Protein.Group": group, "Modified.Sequence": mod,
            "Stripped.Sequence": stripped, "Precursor.Charge": z, "RT": rt, "RT.Start": start,
            "RT.Stop": stop, "Precursor.Quantity": "500"}


class TestDisagreements(unittest.TestCase):
    def test_flags_a_different_peak_only(self):
        sky = sp.skyline_precursors([
            sky_row("sp|P1|A_HUMAN", "PEPC[+57]TIDEK", 2, "r1", 10.00, 9.95, 10.08),   # same peak
            sky_row("sp|P1|A_HUMAN", "PEPC[+57]TIDEK", 2, "r2", 10.20, 10.15, 10.25),  # 12 s away
            sky_row("sp|P1|A_HUMAN", "LMNOPR", 3, "r1", 20.06, 20.00, 20.10),          # 3.6 s, but
            sky_row("sp|P1|A_HUMAN", "LMNOPR", 3, "r2", "#N/A", "#N/A", "#N/A"),       # no peak
            sky_row("sp|P2|B_HUMAN", "PEPC[+57]TIDEK", 2, "r1", 10.00, 9.95, 10.08),   # shared
            sky_row("sp|P1|A_HUMAN", "QQQK", 2, "r1", 5.0, 4.9, 5.1)])                 # no DIA-NN
        dia = sp.diann_index([
            dia_row("P1", "PEPC(UniMod:4)TIDEK", 2, "r1", 10.01, 9.96, 10.07),
            dia_row("P1", "PEPC(UniMod:4)TIDEK", 2, "r2", 10.00, 9.95, 10.06),
            dia_row("P1", "LMNOPR", 3, "r1", 20.00, 19.98, 20.05),
            dia_row("P1", "LMNOPR", 3, "r2", 20.00, 19.95, 20.06)])
        compared, flagged = sp.disagreements(sky, dia)
        self.assertEqual(len(compared), 4)            # the shared peptide counted once
        why = {(r["precursor"], r["run"]): r["reason"] for r in flagged}
        self.assertEqual(set(why), {("PEPCTIDEK/2+", "r2"), ("LMNOPR/3+", "r1"),
                                    ("LMNOPR/3+", "r2")})
        self.assertIn("12.0 s from DIA-NN's", why[("PEPCTIDEK/2+", "r2")])
        self.assertEqual(why[("LMNOPR/3+", "r1")], "apex outside DIA-NN's peak edges")
        self.assertEqual(why[("LMNOPR/3+", "r2")], "Skyline found no peak")
        s = sp.agreement_summary(compared)
        self.assertEqual(s["compared"], 4)

    def test_a_shared_peptide_is_drawn_once(self):
        row = {"FileName": "r1.raw", "PeptideModifiedSequence": "PEPTIDEK", "PrecursorCharge": "2",
               "FragmentIon": "y4", "ProductCharge": "1", "Times": "1,2", "Intensities": "3,4"}
        ch = sp.chromatograms([row, dict(row), dict(row, FragmentIon="precursor")])
        self.assertEqual([lab for lab, _t, _y in ch[(sp.precursor_key("PEPTIDEK", 2), "r1")]],
                         ["y4+"])

    def test_keys_match_across_the_two_spellings(self):
        self.assertEqual(sp.precursor_key("AC[+57]DM[+16]K", 2),
                         sp.precursor_key("AC(UniMod:4)DM(UniMod:35)K", "2"))
        self.assertNotEqual(sp.precursor_key("ACDM[+16]K", 2), sp.precursor_key("ACDMK", 2))
        self.assertEqual(sp.precursor_key("[+42]ACK", 2), sp.precursor_key("(UniMod:1)ACK", 2))
        self.assertEqual(sp.run_of("Z:\\x\\run_1.d"), "run_1")

    def test_table_columns(self):
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "d.tsv")
            sp.write_disagreements(p, [{"protein": "P1", "precursor": "X/2+", "run": "r",
                                        "reason": "Skyline found no peak"}], {"P1": "GENE1"})
            with open(p) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(rows[0]["gene"], "GENE1")
            self.assertEqual(list(rows[0]), sp.DISAGREEMENT_COLUMNS)


OK_REC = {"schema": sp.SCHEMA, "status": "ok", "name": "T",
          "tool": {"name": "Skyline", "version": SKYLINE_PIN["version"]},
          "top": {"value": 10, "source": "skill default"},
          "selection": {"mode": "de_hits", "adjp": 0.05, "adjp_source": "de_provenance.json"},
          "proteins": [{"protein_group": "P1", "gene": "G1", "contrast": "A - B", "adj_p": 1e-4,
                        "drawn": 3, "compared": 12, "flagged": 1}],
          "family": "orbitrap",
          "settings": {"enzyme": "Trypsin/P", "missed_cleavages": 1, "min_length": 7,
                       "max_length": 30, "charges": [2, 3, 4], "fragment_mz": [200, 1800]},
          "mass_accuracy": {"ms1_ppm": 24, "ms2_ppm": 24, "source": "measured by DIA-NN"},
          "runs": {"imported": ["a", "b"], "failed": []},
          "results": {"compared": 12, "flagged": 1, "apex_median_abs_s": 0.8,
                      "worst": [{"gene": "G1", "precursor": "X/2+", "run": "a",
                                 "apex_difference_s": 9.1, "reason": "apex 9.1 s from DIA-NN's"}]},
          "figures": [{"kind": "peaks", "file": "peaks_G1.png", "title": "Peaks of G1 (P1)",
                       "caption": "c"}, {"kind": "areas", "file": "areas_G1.png"}],
          "files": {"precursor_report": "skyline_precursors.tsv",
                    "disagreements": "disagreements.tsv", "share_zip": "T_top_hits.sky.zip"}}


class TestWording(unittest.TestCase):
    """A visual check of DIA-NN's peaks -- never "confirmed by Skyline" or "validated"."""

    def test_the_mandated_words_and_no_overclaim(self):
        sec = sp.report_section(dict(OK_REC))
        md = sec["md"]
        self.assertIn("These pictures are a visual check of DIA-NN's peaks.", md)
        self.assertIn("Skyline reuses DIA-NN's peak edges by default, so a peak that matches is "
                      "not independent confirmation of it.", md)
        self.assertIn("The real warning signs are the peaks where Skyline picked a different "
                      "peak.", md)
        self.assertIn("(DEFAULT — not user-confirmed)", md)     # N = 10 was the skill's default
        self.assertEqual(sec["figures"], ["peaks_G1.png"])
        self.assertIn("![Peaks of G1 (P1)](skyline/peaks_G1.png)", md)
        for status in ("ok", "skipped", "planned"):
            rec = dict(OK_REC, status=status, reason="the container is not reachable here")
            text = sp.report_section(rec)["md"]
            self.assertEqual(client_text.peak_check_overclaims(text), [], status)
            self.assertNotRegex(text.lower(), r"confirmed by|validated|verified by")
        rec = dict(OK_REC)
        rec["methods_paragraph"] = sp.methods_paragraph(rec)
        self.assertEqual(client_text.peak_check_overclaims(rec["methods_paragraph"]), [])
        self.assertIn("not an independent identification", rec["methods_paragraph"])

    def test_the_lint_catches_overclaims(self):
        for bad in ("The peaks were confirmed by Skyline.", "Skyline validated the top hits.",
                    "The peak check verified all ten proteins.",
                    "Skyline independently confirmed the fold change.",
                    "Skyline agreed for 9 proteins, which confirms the quantities."):
            self.assertTrue(client_text.peak_check_overclaims(bad), bad)
        for fine in ("A matching peak is not independent confirmation, as Skyline reuses the "
                     "edges.", "GAPDH was validated by western blot.",
                     "These pictures are a visual check of DIA-NN's peaks."):
            self.assertEqual(client_text.peak_check_overclaims(fine), [], fine)


class TestNotRun(unittest.TestCase):
    """No container: the record says why, and the report, the MANIFEST and the plan's exit do."""

    def make_search(self, root):
        s = os.path.join(root, "search")
        os.makedirs(s)
        raws = []
        for n in ("r1.raw", "r2.raw"):
            raws.append(os.path.join(root, n))
            open(raws[-1], "w").close()
        fasta = os.path.join(root, "db.fasta")
        with open(fasta, "w") as fh:
            fh.write(">sp|P1|A_HUMAN\nPEPTIDEK\n")
        log = AUTO_LOG.replace("/db/search.fasta", fasta)
        for fn, text in (("report.log.txt", log), ("report.parquet", ""),
                         ("report-lib.parquet.skyline.speclib", "")):
            with open(os.path.join(s, fn), "w") as fh:
                fh.write(text)
        with open(os.path.join(s, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "diann", "files": raws, "fasta": fasta}, fh)
        return s

    def test_missing_container(self):
        with tempfile.TemporaryDirectory() as td:
            search = self.make_search(td)
            out = os.path.join(td, "output", "skyline")
            with mock.patch("sys.stdout", new=io.StringIO()) as so:
                rc = sp.main(["plan", "--search", search, "--out", out,
                              "--sif", os.path.join(td, "missing.sif")])
            self.assertEqual(rc, 3)
            self.assertEqual(json.loads(so.getvalue())["status"], "skipped")
            rec = sp.load_record(out)
            self.assertEqual(rec["status"], "skipped")
            self.assertIn("pinned Skyline container is not reachable", rec["reason"])
            sec = sp.report_section(rec)["md"]
            self.assertTrue(sec.startswith("**Not run.**"))
            self.assertIn("not reachable", sec)
            level, name, note = sp.manifest_line(os.path.dirname(out), has_de=True)
            self.assertEqual(level, "SKIPPED")
            self.assertIn("not reachable", note)
            # finalize's step writes exactly that line
            import session
            import make_deposit
            man = make_deposit.Manifest()
            session.peak_check_step({"output_dir": os.path.dirname(out),
                                     "de_dir": os.path.join(td, "output", "tables")}, man)
            self.assertTrue(man.lines[0].startswith("[SKIPPED] Peak check (Skyline"), man.lines)
            self.assertIn("not reachable", man.lines[0])

    def test_no_record_with_a_de_is_skipped_without_a_de_nothing(self):
        with tempfile.TemporaryDirectory() as td:
            self.assertEqual(sp.manifest_line(td, has_de=True)[0], "SKIPPED")
            self.assertIsNone(sp.manifest_line(td, has_de=False))

    def test_turned_off_by_the_user(self):
        with tempfile.TemporaryDirectory() as td, mock.patch("sys.stdout", new=io.StringIO()):
            sp.main(["skip", "--out", td, "--reason", "the lab only wants the tables"])
            rec = sp.load_record(td)
        self.assertIn("Turned off for this analysis: the lab only wants the tables", rec["reason"])

    def test_plan_writes_the_jobs_when_everything_is_there(self):
        with tempfile.TemporaryDirectory() as td:
            search = self.make_search(td)
            sif = os.path.join(td, "sky.sif")
            open(sif, "w").close()
            fake = os.path.join(td, "slurm")
            os.makedirs(fake)
            sb = os.path.join(fake, "sbatch")
            with open(sb, "w") as fh:
                fh.write("#!/bin/sh\nexit 0\n")
            os.chmod(sb, os.stat(sb).st_mode | stat.S_IEXEC)
            out = os.path.join(td, "output", "skyline")
            env = {"HPCCF_SLURM_BIN": fake, "PATH": os.environ.get("PATH", "")}
            with mock.patch.dict(os.environ, env), mock.patch("sys.stdout", new=io.StringIO()) as so, \
                    mock.patch("run_search._associations", return_value=[]):
                os.environ.pop("SLURM_JOB_ID", None)
                rc = sp.main(["plan", "--search", search, "--out", out, "--sif", sif,
                              "--python", sys.executable, "--partition", "high", "--account",
                              "acct", "--qos", "acct-high-qos", "--name", "T"])
            self.assertEqual(rc, 0, so.getvalue())
            w = os.path.join(out, sp.WORK)
            with open(os.path.join(w, "sky_2_import.sbatch")) as fh:
                arr = fh.read()
            self.assertIn("#SBATCH --array=1-2%2", arr)          # one task per raw file
            with open(os.path.join(w, "plan.json")) as fh:
                plan = json.load(fh)
            self.assertEqual(plan["mass_accuracy"]["ms2_ppm"], 24.0)
            self.assertEqual(plan["top"], {"value": sp.TOP_N_DEFAULT, "source": "skill default"})
            self.assertEqual(sp.load_record(out)["status"], "planned")
            self.assertTrue(os.path.islink(os.path.join(w, "search",
                                                        "report-lib.parquet.skyline.speclib")))
            self.assertIn("submit with", json.loads(so.getvalue()))
            # a second plan before the first has a result is refused; --replan clears the old
            # run's outputs, caches and status files
            with open(os.path.join(w, "status", "import_1.json"), "w") as fh:
                fh.write('{"rc": 0}')
            with open(os.path.join(out, "peaks_OLD.png"), "w") as fh:
                fh.write("x")
            with mock.patch.dict(os.environ, env), mock.patch("sys.stdout", new=io.StringIO()), \
                    mock.patch("run_search._associations", return_value=[]):
                with self.assertRaises(SystemExit):
                    sp.main(["plan", "--search", search, "--out", out, "--sif", sif,
                             "--python", sys.executable, "--partition", "high", "--account",
                             "acct", "--qos", "acct-high-qos"])
                rc = sp.main(["plan", "--search", search, "--out", out, "--sif", sif,
                              "--python", sys.executable, "--partition", "high", "--account",
                              "acct", "--qos", "acct-high-qos", "--replan"])
            self.assertEqual(rc, 0)
            self.assertFalse(os.path.exists(os.path.join(w, "status", "import_1.json")))
            self.assertFalse(os.path.exists(os.path.join(out, "peaks_OLD.png")))


class TestBrief(unittest.TestCase):
    """The analysis brief tells the report's writer the page adds the section, and the words."""

    def test_brief(self):
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            tables, figs, sky = (os.path.join(td, x) for x in ("tables", "figures", "skyline"))
            for d in (tables, figs, sky):
                os.makedirs(d)
            with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
                json.dump({"adjp": 0.05}, fh)
            with open(os.path.join(sky, sp.RECORD), "w") as fh:
                json.dump(sp.skipped_record("the container is not reachable here", "plan"), fh)
            out = os.path.join(td, "brief.md")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"),
                                "--out", out, "--de-dir", tables, "--figures-dir", figs],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(out, encoding="utf-8") as fh:
                brief = fh.read()
        self.assertIn('Do not write a "Peak check (Skyline)" section', brief)
        self.assertIn("it was not run: the container is not reachable here", brief)
        self.assertIn("a matching peak is not independent confirmation", brief)


class TestFollowUp(unittest.TestCase):
    """PRM / SRM lists: scheduled only on the same instrument family AND gradient."""

    REC = {"family": "orbitrap", "status": "ok"}

    def test_scheduling_rules(self):
        t = {x["key"]: x for x in sp.prm_plan(self.REC, ["thermo_fusion", "bruker_timstof"],
                                               True, 60)}
        self.assertTrue(t["thermo_fusion"]["scheduled"])
        self.assertIn("--exp-method-type=Scheduled", t["thermo_fusion"]["extra"])
        self.assertFalse(t["bruker_timstof"]["scheduled"])            # another family
        self.assertIn("do not transfer", t["bruker_timstof"]["why"])
        self.assertIn("--exp-method-type=Standard", t["bruker_timstof"]["extra"])
        t = sp.prm_plan(self.REC, ["thermo_qe"], False, 45)[0]       # gradient not confirmed
        self.assertFalse(t["scheduled"])
        self.assertIn("--exp-run-length=45", t["extra"])
        with self.assertRaises(ValueError):                           # Thermo, unscheduled
            sp.prm_plan(self.REC, ["thermo_qe"], False, None)
        with self.assertRaises(ValueError) as e:
            sp.prm_plan(self.REC, ["thermo_astral"], True, 60)
        self.assertIn("no isolation list", str(e.exception))
        with self.assertRaises(ValueError):
            sp.prm_plan(self.REC, ["sciex"], True, 60)

    def test_the_batch(self):
        p = plan_for()
        targets = sp.prm_plan({"family": "timstof"}, ["bruker_timstof", "thermo_srm"], True, 30)
        b = sp.prm_batch(p, p["doc"], targets, "/out/followup")
        lines = b.splitlines()
        self.assertIn("--refine-max-peptide-peak-rank=3", lines[0])
        self.assertIn("--refine-max-transition-peak-rank=6", lines[0])
        self.assertIn("--refine-max-precursor-only", lines[0])
        self.assertIn("--full-scan-acquisition-method=PRM", lines[1])
        self.assertIn("--full-scan-precursor-isotopes=None", lines[1])
        self.assertIn('--exp-isolationlist-instrument="Bruker timsTOF"', lines[2])
        self.assertIn("--exp-method-type=Scheduled", lines[2])
        srm = [x for x in lines if "--tran-predict-ce" in x][0]
        self.assertIn('--tran-product-ion-types="y, b"', srm)          # the "p" ion type dropped
        self.assertIn('--tran-predict-ce="Thermo TSQ Quantiva"', srm)
        self.assertIn('--exp-translist-instrument="Thermo Quantiva"', b)
        self.assertNotIn("import-process-count", b)
        for _, name, _ in sp.PRM_INSTRUMENTS.values():      # only what SkylineCmd 26.1 offers
            self.assertIn(name, ("Bruker timsTOF", "Thermo Fusion", "Thermo Q Exactive",
                                 "Thermo Stellar", "Thermo Quantiva"))

    def test_readme_says_ce_is_left_for_the_method(self):
        rec = {"name": "T", "family": "timstof"}
        spec = {"targets": sp.prm_plan(rec, ["bruker_timstof"], True, None)}
        text = sp.followup_readme(spec, rec, {spec["targets"][0]["file"]: 15})
        self.assertIn("leaves the collision energy at 0", text)
        self.assertIn("1/K0 window measured", text)
        rec = {"name": "T", "family": "orbitrap"}           # no ion mobility in the search
        spec = {"targets": sp.prm_plan(rec, ["bruker_timstof"], True, None)}
        text = sp.followup_readme(spec, rec, {spec["targets"][0]["file"]: 15})
        self.assertIn("Skyline's default, NOT a measurement", text)


class TestPinsAndPlumbing(unittest.TestCase):
    def test_the_pin(self):
        self.assertEqual(SKYLINE_PIN["version"], "26.1.0.266")
        self.assertEqual(SKYLINE_PIN["sif_sha256"],
                         "b5f135587e9b589dd2d4686ff520573ea237b84f53b1343d5b76213c16f302bf")
        self.assertTrue(SKYLINE_PIN["sif"].endswith(SKYLINE_PIN["tag"] + ".sif"))

    def test_container_sha_is_checked(self):
        with tempfile.TemporaryDirectory() as td:
            fake = os.path.join(td, os.path.basename(SKYLINE_PIN["sif"]))
            with open(fake, "w") as fh:
                fh.write("not the image")
            p = dict(plan_for(work=os.path.join(td, "skyline_work")), sif=fake)
            with self.assertRaises(sp.Skip) as e:
                sp.verify_container(p)
            self.assertIn("does not match its pinned sha256", str(e.exception))

    def test_work_folder_is_scratch(self):
        self.assertTrue(scratch_files.is_scratch_dir(sp.WORK))

    def test_output_files_catalog(self):
        import make_report
        for n in ("peaks_GAPDH.png", "areas_GAPDH.png", "disagreements.tsv",
                  "skyline_precursors.tsv", "X_top_hits.sky.zip", "skyline_provenance.json"):
            self.assertEqual(make_report.describe(n)[0], make_report.PEAKS, n)

    def test_deliver_takes_the_folder_not_the_work(self):
        import core_submission as cs
        with tempfile.TemporaryDirectory() as td:
            os.makedirs(os.path.join(td, "skyline", sp.WORK))
            for f in ("skyline/peaks_G.png", "skyline/skyline_work/doc.skyd",
                      "Analysis_Report.html"):
                open(os.path.join(td, f), "w").close()
            names = {i["name"]: i["status"] for i in cs.plan_delivery(td)}
        self.assertEqual(names.get("skyline/peaks_G.png"), "OK")
        self.assertEqual(names.get("skyline/skyline_work/"), "SKIPPED")
        self.assertNotIn("skyline/skyline_work/doc.skyd", names)


@unittest.skipUnless(HAVE_MPL, "needs matplotlib")
class TestPictures(unittest.TestCase):
    def test_finish_draws_and_records(self):
        with tempfile.TemporaryDirectory() as td:
            out = os.path.join(td, "skyline")
            w = os.path.join(out, sp.WORK)
            os.makedirs(w)
            os.environ.setdefault("MPLCONFIGDIR", os.path.join(td, ".mpl"))
            sky = [sky_row("sp|P1|A_HUMAN", "PEPTIDEK", 2, r, 10.0 + d, 9.95 + d, 10.06 + d)
                   for r, d in (("r1", 0), ("r2", 0.2))]
            with open(os.path.join(out, "skyline_precursors.tsv"), "w", newline="") as fh:
                wr = csv.DictWriter(fh, SKY_HEADER, delimiter="\t")
                wr.writeheader()
                wr.writerows(sky)
            with open(os.path.join(w, "diann_targets.tsv"), "w", newline="") as fh:
                wr = csv.DictWriter(fh, ["Run", "Protein.Group", "Modified.Sequence",
                                         "Stripped.Sequence", "Precursor.Charge", "RT",
                                         "RT.Start", "RT.Stop", "Precursor.Quantity"],
                                    delimiter="\t")
                wr.writeheader()
                for r in ("r1", "r2"):
                    wr.writerow(dia_row("P1", "PEPTIDEK", 2, r, 10.0, 9.95, 10.06))
            times = ",".join(f"{9.8 + i * 0.02:.2f}" for i in range(20))
            ints = ",".join(str(i * (20 - i)) for i in range(20))
            with open(os.path.join(w, "skyline_chromatograms.tsv"), "w", newline="") as fh:
                wr = csv.writer(fh, delimiter="\t")
                wr.writerow(["FileName", "PeptideModifiedSequence", "PrecursorCharge",
                             "ProductMz", "FragmentIon", "ProductCharge", "Times", "Intensities"])
                for r in ("r1", "r2"):
                    for ion in ("y4", "y5", "precursor"):
                        wr.writerow([r + ".raw", "PEPTIDEK", 2, 500, ion, 1, times, ints])
            with open(os.path.join(w, "skyline_transitions.tsv"), "w", newline="") as fh:
                wr = csv.writer(fh, delimiter="\t")
                wr.writerow(["PeptideModifiedSequence", "PrecursorCharge", "FragmentIon",
                             "ProductCharge", "FileName", "Area"])
                for r in ("r1", "r2"):
                    wr.writerow(["PEPTIDEK", 2, "y4", 1, r + ".raw", 100])
            p = dict(plan_for(raw=["/d/r1.raw", "/d/r2.raw"], work=w), out=out,
                     share_zip=os.path.join(out, "T_top_hits.sky.zip"), group_of={},
                     family_source="test", top={"value": 1, "source": "user"},
                     search={"dir": "/s"}, skill_version="t")
            build = {"tool": {"name": "Skyline", "version": SKYLINE_PIN["version"]},
                     "selection": {"mode": "de_hits", "adjp": 0.05,
                                   "adjp_source": "de_provenance.json",
                                   "targets": [{"protein_group": "P1", "gene": "GENE1",
                                                "contrast": "A - B", "adj_p": 0.001}]},
                     "rt_window": {"minutes": 1.0, "source": "t"}}
            rec = sp.finish(p, build, p["raw_files"], [])
            self.assertEqual(rec["results"]["flagged"], 1)               # r2: 12 s away
            self.assertTrue(os.path.getsize(os.path.join(out, "peaks_GENE1.png")) > 1000)
            self.assertTrue(os.path.getsize(os.path.join(out, "areas_GENE1.png")) > 1000)
            self.assertTrue(os.path.isfile(os.path.join(out, "peaks_disagreements.png")))
            self.assertEqual(client_text.peak_check_overclaims(rec["methods_paragraph"]), [])


if __name__ == "__main__":
    unittest.main()
