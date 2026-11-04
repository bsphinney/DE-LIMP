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
import getpass
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
sys.path.insert(0, HERE)

import client_text      # noqa: E402
import scratch_files    # noqa: E402
import skyline_peaks as sp  # noqa: E402
from resolve_defaults import SKYLINE_PIN  # noqa: E402

try:
    import matplotlib  # noqa: F401  -- only whether it is there
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False
try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:                       # the stdlib-only CI: those tests are skipped
    pa = pq = None

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
        self.assertTrue(any("precursor m/z range (357-1105)" in x for x in s["not_carried_over"]))

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
        # no fixed cysteine modification is a setting the search states, not a missing one
        # (2.14.2, I37): the peak check runs, with Skyline's default carbamidomethyl cleared
        self.assertEqual(sp.search_settings(sp.command_flags(CMD.replace("--unimod4 ", "")))
                         ["fixed_mods"], [])

    def test_no_cam_search_clears_skyline_default_mods(self):
        """2.14.2, I37 (I25): unalkylated samples are searched with no fixed C modification, but
        SkylineCmd --new starts with Skyline's default Carbamidomethyl (C), so no cysteine peptide
        of DIA-NN's library would match. --pep-clear-mods (SkylineCmd 26.1.0.266 --help, HIVE job
        25386172) clears it, right after --new and --overwrite."""
        s = sp.search_settings(sp.command_flags(CMD.replace("--unimod4 ", "")))
        self.assertEqual(s["fixed_mods"], [])
        self.assertEqual(s["mod_args"], ["--pep-clear-mods"])
        p = dict(plan_for(), settings=s)
        args = sp.build_args(p, 2.0)
        self.assertTrue(args[0].startswith("--new="))
        self.assertEqual(args[1:3], ["--overwrite", "--pep-clear-mods"])
        self.assertEqual(args.count("--pep-clear-mods"), 1)
        # recorded: the record carries the search's settings, and what was done to Skyline's
        p.update(plan_id="p", skill_version="t", top={"value": 10, "source": "skill default"},
                 family_source="t", search={"dir": "/s"})
        self.assertEqual(sp.planned_record(p)["settings"]["mod_args"], ["--pep-clear-mods"])
        rec = dict(OK_REC, settings=dict(OK_REC["settings"], fixed_mods=[]))
        self.assertIn("fragments 200-1800 m/z, no fixed modification, as searched: cysteine "
                      "unmodified, and Skyline's default carbamidomethyl modification cleared), ",
                      sp.methods_paragraph(rec))

    def test_methods_name_the_cam_fixed_modification(self):
        """(pin, L9 review F8) every 2.14.2 --unimod4 record carries fixed_mods: its Methods name
        carbamidomethyl as the fixed modification and say nothing was cleared."""
        s = sp.search_settings(sp.command_flags(CMD))
        rec = dict(OK_REC, settings=dict(OK_REC["settings"], fixed_mods=s["fixed_mods"]))
        text = sp.methods_paragraph(rec)
        self.assertIn("fragments 200-1800 m/z, Carbamidomethyl (C) as a fixed modification), "
                      "centroided spectra", text)
        self.assertNotIn("cleared", text)
        self.assertNotIn("unmodified", text)

    def test_methods_with_a_variable_cysteine_modification(self):
        """L9 review F8: a search with a variable C modification and no fixed one (--var-mod
        UniMod:4,...,C without --unimod4) also gets Skyline's default cleared, but its cysteine
        is not "unmodified". And no variant closes two parentheses at once ("))")."""
        s = sp.search_settings(sp.command_flags(CMD.replace(
            "--unimod4 ", "--var-mod UniMod:4,57.021464,C --var-mod UniMod:35,15.994915,M ")))
        self.assertEqual((s["fixed_mods"], s["mod_args"]), ([], ["--pep-clear-mods"]))
        rec = dict(OK_REC, settings=dict(OK_REC["settings"], fixed_mods=[],
                                         var_mods=s["var_mods"]))
        text = sp.methods_paragraph(rec)
        self.assertIn("fragments 200-1800 m/z, no fixed modification, as searched: Skyline's "
                      "default carbamidomethyl modification cleared), ", text)
        self.assertNotIn("unmodified", text)
        # only a variable modification ON C counts: Met-ox alone leaves cysteine unmodified
        rec["settings"]["var_mods"] = ["UniMod:35,15.994915,M"]
        self.assertIn("cysteine unmodified", sp.methods_paragraph(rec))
        for fixed in ([], ["Carbamidomethyl (C)"]):
            for var in ([], ["UniMod:4,57.021464,C"]):
                r = dict(OK_REC, settings=dict(OK_REC["settings"], fixed_mods=fixed, var_mods=var))
                self.assertNotIn("))", sp.methods_paragraph(r), (fixed, var))

    def test_methods_name_the_variable_mods_and_what_was_not_carried_over(self):
        """H4 finding 3 (Ruling 76): the Methods paragraph left out the search's variable
        modifications and the settings not carried over (N-terminal Met excision, the precursor
        m/z range), which the record holds. It names them from the record's settings; a search
        with none says so; a record without the field says nothing (not guessed)."""
        s = sp.search_settings(sp.command_flags(CMD + " --var-mod UniMod:35,15.994915,M "
                                                      "--var-mod UniMod:1,42.010565,*n --var-mods 2"))
        self.assertEqual(len(s["not_carried_over"]), 2)            # --met-excision, --min/max-pr-mz
        text = sp.methods_paragraph(dict(OK_REC, settings=dict(OK_REC["settings"], **{
            k: s[k] for k in ("var_mods", "max_var_mods", "not_carried_over", "fixed_mods")})))
        self.assertIn("The modifications found in DIA-NN's library were added to the document; "
                      "the search's variable modifications were Oxidation (M); Acetyl (Protein "
                      "N-term), at most 2 per peptide. ", text)
        self.assertIn("Not carried over from the search: DIA-NN's precursor m/z range "
                      "(357-1105)", text)
        self.assertIn("N-terminal methionine excision (--met-excision)", text)
        self.assertNotIn("))", text)
        # the same words as the report's (one definition)
        nco = sp.not_carried_over_sentence(s)
        self.assertIn(nco, text)
        self.assertIn(nco, sp.report_section(dict(OK_REC, settings=dict(OK_REC["settings"], **s)))["md"])
        # no variable modification: said; no field (a record that does not say): nothing
        none = sp.methods_paragraph(dict(OK_REC, settings=dict(OK_REC["settings"], var_mods=[])))
        self.assertIn("The search had no variable modifications.", none)
        self.assertNotIn("Not carried over", none)
        old = sp.methods_paragraph(dict(OK_REC))
        self.assertNotIn("variable modification", old)
        self.assertNotIn("Not carried over", old)

    def test_other_fixed_c_mod_still_skips(self):
        for flags in ("--fixed-mod UniMod:24,71.037114,C ", "--unimod4 --fixed-mod "
                      "UniMod:39,45.987721,C ", "--fixed-mod Carbamidomethyl,57.021464,C "):
            with self.subTest(flags=flags), self.assertRaises(sp.Skip) as e:
                sp.search_settings(sp.command_flags(CMD.replace("--unimod4 ", flags)))
            self.assertIn("a fixed cysteine modification other than carbamidomethyl is not yet "
                          "tested with Skyline", str(e.exception))
            self.assertIn(flags.split("--fixed-mod ")[1].strip(), str(e.exception))
        # a fixed modification on another residue is not a cysteine one
        s = sp.search_settings(sp.command_flags(CMD + " --fixed-mod SILAC,0.0,KR,label"))
        self.assertEqual(s["fixed_mods"], ["Carbamidomethyl (C)"])

    def test_cam_search_unchanged(self):
        """(pin) --unimod4: Skyline's default Carbamidomethyl (C) stays; nothing is cleared."""
        s = sp.search_settings(sp.command_flags(CMD))
        self.assertEqual(s["fixed_mods"], ["Carbamidomethyl (C)"])
        args = sp.build_args(dict(plan_for(), settings=s), 2.0)
        self.assertNotIn("--pep-clear-mods", args)
        self.assertEqual(args[1], "--overwrite")

    def test_cam_as_a_fixed_mod_is_cam(self):
        s = sp.search_settings(sp.command_flags(CMD.replace(
            "--unimod4 ", "--fixed-mod UniMod:4,57.021464,C ")))
        self.assertEqual(s["fixed_mods"], ["Carbamidomethyl (C)"])
        self.assertNotIn("--pep-clear-mods", sp.build_args(dict(plan_for(), settings=s), 2.0))

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

            def fake(plan, args, log_name, rw=(), timeout_s=None):
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
            sp.run_skyline(p, ["--import-process-count=2"], "x.txt")

    def test_join_batch_lists_only_the_imported_files(self):
        p = plan_for()
        batch = sp.join_batch(p, p["raw_files"][:2])
        first = batch.splitlines()[0]
        self.assertEqual(first.count("--import-file="), 2)
        self.assertNotIn("r3.raw", batch)
        self.assertIn('--in="Z:\\w\\skyline_work\\T_top_hits.sky"', first)
        self.assertIn("--share-type=minimal", batch)
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
        self.assertEqual(set(why), {("PEPC[+57]TIDEK/2+", "r2"), ("LMNOPR/3+", "r1"),
                                    ("LMNOPR/3+", "r2")})
        self.assertIn("12.0 s from DIA-NN's", why[("PEPC[+57]TIDEK/2+", "r2")])
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
            self.assertEqual(client_text.peak_check_overclaims_loose(text), [], status)
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
        # review of 2.14 (LOW-2): the verb is another check's, or the reader is asked to do it --
        # loose (a warning for a person), never strict
        for honest in ("The peak check drew the actin peaks, and the contaminant check confirmed "
                       "that bovine actin took ACTB peptides.",
                       "Open the Skyline file to verify the peaks yourself.",
                       "Looking at the Skyline pictures, we confirmed the peak edges by eye."):
            with self.subTest(sentence=honest):
                self.assertEqual(client_text.peak_check_overclaims(honest), [])
                self.assertTrue(client_text.peak_check_overclaims_loose(honest))
        for bad in ("The peaks were confirmed by Skyline.",
                    "Skyline agreed for 9 proteins, which confirms the quantities.",
                    "The contaminant check and Skyline both confirmed the ACTB loss."):
            with self.subTest(sentence=bad):
                self.assertTrue(client_text.peak_check_overclaims(bad))
                self.assertEqual(client_text.peak_check_overclaims_loose(bad), [])


@unittest.skipUnless(pq, "pyarrow is not installed")
class TestQFilter(unittest.TestCase):
    """Review of 2.14 (MED-1): the peak check had its own copy of run_de.R's q filter, and it
    dropped a NaN or missing q-value that the DE keeps. It now reads through
    diann_q_columns.keep_expression, the one definition."""

    def test_nan_and_null_stay_as_in_the_de(self):
        with tempfile.TemporaryDirectory() as td:
            rep = os.path.join(td, "report.parquet")
            ok = [0.001] * 4
            pq.write_table(pa.table({
                "Run": ["r1"] * 4, "Protein.Group": ["P1"] * 4,
                "Precursor.Id": ["a", "b", "c", "d"],
                "Q.Value": pa.array([0.001, float("nan"), None, 0.5], pa.float64()),
                "Lib.Q.Value": ok, "Lib.PG.Q.Value": ok, "PG.Q.Value": ok}), rep)
            t, cols = sp._q_filtered(rep, ["Run", "Precursor.Id", "Missing.Column"], 0.01)
        self.assertEqual(cols, ["Run", "Precursor.Id"])
        self.assertEqual(sorted(t["Precursor.Id"].to_pylist()), ["a", "b", "c"])

    def test_pg_q_value_keeps_its_own_cutoff(self):
        with tempfile.TemporaryDirectory() as td:
            rep = os.path.join(td, "report.parquet")
            pq.write_table(pa.table({"Precursor.Id": ["a", "b"], "Q.Value": [0.001, 0.001],
                                     "PG.Q.Value": [0.04, 0.06]}), rep)
            t, _cols = sp._q_filtered(rep, ["Precursor.Id"], 0.01)
        self.assertEqual(t["Precursor.Id"].to_pylist(), ["a"])     # 0.05, not 0.01


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
            got = json.loads(so.getvalue())
            self.assertIn("submit with", got)
            # review of 2.14 (LOW-1): the board's job record carries the figure a script printed
            # -- build + import (one task per raw file) + join, CPUs x time limit, summed
            import cpu_hours
            res = cpu_hours.reservation([os.path.join(w, f"sky_{i}.sbatch") for i in
                                         ("1_build", "2_import", "3_join")])
            self.assertEqual(got["cpu_hours_reserved"], res["cpu_hours_reserved"])
            self.assertGreater(got["cpu_hours_reserved"], 0)
            self.assertEqual(got["gpu_hours_reserved"], 0)
            with open(os.path.join(w, "submit.sh")) as fh:
                lines = fh.read().splitlines()
            self.assertEqual(lines[0], "#!/bin/bash")
            self.assertEqual(lines[1], cpu_hours.submit_comment(res))
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


def plan_here(td, search, *extra):
    """`plan` with everything there (a stand-in sbatch, a container file) -> (rc, out, stdout)."""
    sif = os.path.join(td, "sky.sif")
    open(sif, "w").close()
    fake = os.path.join(td, "slurm")
    os.makedirs(fake, exist_ok=True)
    sb = os.path.join(fake, "sbatch")
    with open(sb, "w") as fh:
        fh.write("#!/bin/sh\nexit 0\n")
    os.chmod(sb, os.stat(sb).st_mode | stat.S_IEXEC)
    out = os.path.join(td, "output", "skyline")
    env = {"HPCCF_SLURM_BIN": fake, "PATH": os.environ.get("PATH", "")}
    with mock.patch.dict(os.environ, env), mock.patch("sys.stdout", new=io.StringIO()) as so, \
            mock.patch("run_search._associations", return_value=[]):
        os.environ.pop("SLURM_JOB_ID", None)
        rc = sp.main(["plan", "--search", search, "--out", out, "--sif", sif, "--python",
                      sys.executable, "--partition", "high", "--account", "acct", "--qos",
                      "acct-high-qos", "--name", "T"] + list(extra))
    return rc, out, so.getvalue()


class ChainLayout(unittest.TestCase):
    """2.14.2, I38 (I19): run_search's 5-step chain -- the route for more than 5 DIA-NN files on
    SLURM -- writes its library as empirical.parquet.skyline.speclib (step 3's --out-lib
    empirical.parquet). Skyline pairs a library with the report by a shared leading name ("unable
    to determine DIA-NN report filename for 'empirical.parquet.skyline.speclib'"), so the build
    failed on every chain search. The library is linked under the report's paired name whatever it
    is called itself; a link under that name built and imported on HIVE (jobs 25168317-9)."""

    def test_the_paired_name(self):
        self.assertEqual(sp.speclib_link_name("/s/report.parquet"),
                         "report-lib.parquet.skyline.speclib")
        self.assertEqual(sp.speclib_link_name("/s/x.y.parquet"), "x.y-lib.parquet.skyline.speclib")

    def test_chain_library_is_linked_under_the_report_paired_name(self):
        with tempfile.TemporaryDirectory() as td:
            search = TestNotRun().make_search(td)
            empirical = os.path.join(search, "empirical.parquet.skyline.speclib")
            os.rename(os.path.join(search, "report-lib.parquet.skyline.speclib"), empirical)
            self.assertEqual(sorted(f for f in os.listdir(search) if "parquet" in f),
                             ["empirical.parquet.skyline.speclib", "report.parquet"])
            rc, out, said = plan_here(td, search)
            self.assertEqual(rc, 0, said)
            w = os.path.join(out, sp.WORK)
            with open(os.path.join(w, "plan.json")) as fh:
                plan = json.load(fh)
            link = plan["speclib_link"]
            self.assertEqual(link, os.path.join(w, "search", "report-lib.parquet.skyline.speclib"))
            self.assertTrue(os.path.islink(link))
            self.assertEqual(os.path.realpath(link), os.path.realpath(empirical))
            self.assertEqual(plan["search"]["speclib"], empirical)     # the record names the file
            self.assertFalse(os.path.lexists(os.path.join(w, "search", os.path.basename(empirical))))
            self.assertIn(f"--import-search-file={sp.win(link)}", sp.build_args(plan, 1.0))


class ReusedSearchFolder(unittest.TestCase):
    """L9 review F11: a folder searched single-shot and then re-searched as a chain holds both
    report-lib.parquet.skyline.speclib (the first search's; the chain's clear_stale leaves it)
    and the chain's empirical.parquet.skyline.speclib. The paired name alone took the earlier
    search's library. The search's record (search_provenance.json search_mode) now says which
    library is this search's own, and a choice it cannot settle is refused, naming both."""

    PAIRED = "report-lib.parquet.skyline.speclib"
    CHAIN = "empirical.parquet.skyline.speclib"

    def folder(self, td, files):
        s = os.path.join(td, "s")
        os.makedirs(s)
        for f in files:
            open(os.path.join(s, f), "w").close()
        return s, os.path.join(s, "report.parquet")

    def test_the_chain_name_is_the_chains_own(self):
        import diann_parallel
        self.assertEqual(diann_parallel.EMPIRICAL_LIBRARY + sp.SKYLINE_SPECLIB, self.CHAIN)

    def test_a_chain_search_takes_its_own_library_in_a_reused_folder(self):
        with tempfile.TemporaryDirectory() as td:
            search = TestNotRun().make_search(td)          # holds report-lib...speclib
            empirical = os.path.join(search, self.CHAIN)
            open(empirical, "w").close()
            prov = os.path.join(search, "search_provenance.json")
            with open(prov) as fh:
                rec = json.load(fh)
            with open(prov, "w") as fh:
                json.dump(dict(rec, search_mode="parallel_5step"), fh)
            rc, out, said = plan_here(td, search)
            self.assertEqual(rc, 0, said)
            with open(os.path.join(out, sp.WORK, "plan.json")) as fh:
                plan = json.load(fh)
        self.assertEqual(plan["search"]["speclib"], empirical)
        self.assertIn("search_mode parallel_5step", plan["search"]["speclib_chosen_by"])
        self.assertEqual(os.path.basename(plan["speclib_link"]), self.PAIRED)   # T27's pairing

    def test_each_route_takes_its_own(self):
        with tempfile.TemporaryDirectory() as td:
            s, rep = self.folder(td, [self.PAIRED, self.CHAIN])
            self.assertEqual(sp.find_speclib(s, rep, "parallel_5step")[0],
                             os.path.join(s, self.CHAIN))
            self.assertEqual(sp.find_speclib(s, rep, "single_shot")[0],
                             os.path.join(s, self.PAIRED))

    def test_an_unsettled_choice_is_refused_naming_both(self):
        with tempfile.TemporaryDirectory() as td:
            s, rep = self.folder(td, [self.PAIRED, self.CHAIN])
            with self.assertRaises(sp.Skip) as e:
                sp.find_speclib(s, rep, None)                # no search_mode recorded
        for f in (self.PAIRED, self.CHAIN, "no search_mode"):
            self.assertIn(f, str(e.exception))

    def test_a_route_whose_library_is_gone_never_takes_the_other(self):
        for mode, there, gone in (("parallel_5step", self.PAIRED, self.CHAIN),
                                  ("single_shot", self.CHAIN, self.PAIRED)):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as td:
                s, rep = self.folder(td, [there])
                with self.assertRaises(sp.Skip) as e:
                    sp.find_speclib(s, rep, mode)
                for f in (there, gone, mode):
                    self.assertIn(f, str(e.exception))

    def test_one_library_and_no_record_is_taken_as_before(self):
        # T27 (I38): a chain's folder with no search_mode, and a single search's
        for files in ([self.CHAIN], [self.PAIRED], [self.PAIRED, "x.predicted.speclib"
                                                     + sp.SKYLINE_SPECLIB]):
            with self.subTest(files=files), tempfile.TemporaryDirectory() as td:
                s, rep = self.folder(td, files)
                self.assertEqual(sp.find_speclib(s, rep)[0], os.path.join(s, files[0]))


class TestReviewFindings(unittest.TestCase):
    """What a review of the branch found, each held here."""

    def test_dda_is_not_planned_as_dia(self):
        with tempfile.TemporaryDirectory() as td:
            search = TestNotRun().make_search(td)
            log = os.path.join(search, "report.log.txt")
            with open(log) as fh:
                text = fh.read()
            with open(log, "w") as fh:
                fh.write(text.replace("--unimod4", "--unimod4 --dda"))
            out = os.path.join(td, "output", "skyline")
            with mock.patch("sys.stdout", new=io.StringIO()):
                self.assertEqual(sp.main(["plan", "--search", search, "--out", out]), 3)
            self.assertIn("DDA", sp.load_record(out)["reason"])

    def test_output_folder_inside_a_raw_folder_stays_writable(self):
        # the session's default place is inside the raw-data folder: its read-only bind must come
        # first, or it covers the read-write output folder below it
        p = dict(plan_for(), out="/raw/proj/session/output/skyline",
                 ro_binds=["/raw/proj", "/search"])
        argv = sp.apptainer_argv(p)
        binds = [argv[i + 1] for i, x in enumerate(argv) if x == "--bind"]
        self.assertLess(binds.index("/raw/proj:/raw/proj:ro"),
                        binds.index("/raw/proj/session/output/skyline:"
                                    "/raw/proj/session/output/skyline"))

    def test_a_raw_file_searched_through_a_link(self):
        # DIA-NN cuts its command at '--', so such a run was searched through a link, whose name
        # is the run's name; Skyline imports the link, and both folders are bound
        prov = {"files": ["/raw/a--b.d", "/raw/c.d"],
                "input_links": [{"original": "/raw/a--b.d", "link": "/s/input_links/a_b.d",
                                 "run_name": "a_b", "renamed": True}]}
        raw = sp.raw_files_of(prov, None)
        self.assertEqual(raw, ["/s/input_links/a_b.d", "/raw/c.d"])
        with tempfile.TemporaryDirectory() as td:
            real = os.path.join(td, "real")
            os.makedirs(real)
            open(os.path.join(real, "x.raw"), "w").close()
            os.makedirs(os.path.join(td, "links"))
            os.symlink(os.path.join(real, "x.raw"), os.path.join(td, "links", "x.raw"))
            binds = sp.ro_binds(os.path.join(td, "s"), [os.path.join(td, "links", "x.raw")])
            rt = os.path.realpath(td)
            self.assertTrue(any(b in binds for b in (os.path.join(td, "links"),
                                                      os.path.join(rt, "links"))), binds)
            self.assertIn(os.path.join(rt, "real"), binds)

    def test_two_mods_on_one_residue(self):
        k = sp.precursor_key
        self.assertEqual(k("M[+58]ACK", 2), k("(UniMod:1)M(UniMod:35)ACK", 2))
        self.assertEqual(k("M[+42]ACK", 2), k("(UniMod:1)MACK", 2))
        # one of the two is never the residue with both (an acetyl-only peptide is not its
        # acetyl + oxidation form)
        self.assertNotEqual(k("(UniMod:1)MACK", 2), k("(UniMod:1)M(UniMod:35)ACK", 2))
        self.assertEqual(k("PEPM[+16]K", 2), k("PEPM(UniMod:35)K", 2))
        self.assertEqual(k("PEPM[+15.994915]K", 2), k("PEPM(UniMod:35)K", 2))
        self.assertNotEqual(k("PEPM[+16]K", 2), k("PEPMK", 2))
        self.assertNotEqual(k("PEPM(UniMod:99999)K", 2), k("PEPM[+16]K", 2))     # unknown id

    def write_record(self, out_dir, rec):
        os.makedirs(os.path.join(out_dir, "skyline"), exist_ok=True)
        with open(os.path.join(out_dir, "skyline", sp.RECORD), "w") as fh:
            json.dump(rec, fh) if isinstance(rec, dict) else fh.write(rec)

    def test_a_de_run_again_makes_the_record_out_of_date(self):
        with tempfile.TemporaryDirectory() as td:
            os.makedirs(os.path.join(td, "tables"))
            with open(os.path.join(td, "tables", "de_provenance.json"), "w") as fh:
                fh.write('{"adjp": 0.05}')
            de = os.path.join(td, "tables", "de_provenance.json")
            self.write_record(td, dict(OK_REC, de_record=de, de_record_sha256=sp.de_fingerprint(de)))
            self.assertEqual(sp.manifest_line(td, True)[0], "OK")
            with open(os.path.join(td, "tables", "de_provenance.json"), "w") as fh:
                fh.write('{"adjp": 0.05, "rerun": true}')
            self.assertEqual(sp.record_state(td)[0], "stale")
            # the client's words: what the report holds, never a task for staff (C2 M5)
            self.assertEqual(sp.report_data(td)["md"], sp.REPORT_WORDS["stale"])
            self.assertNotIn("run again", sp.report_data(td)["md"])
            level, _n, note = sp.manifest_line(td, True)
            self.assertEqual(level, "SKIPPED")
            self.assertIn("run again", note)

    def test_a_check_made_without_a_de_goes_stale_when_one_appears(self):
        with tempfile.TemporaryDirectory() as td:
            self.write_record(td, dict(OK_REC, de_record=None, de_record_sha256=None))
            self.assertEqual(sp.record_state(td)[0], "ok")
            os.makedirs(os.path.join(td, "tables"))
            with open(os.path.join(td, "tables", "de_provenance.json"), "w") as fh:
                fh.write("{}")
            self.assertEqual(sp.record_state(td)[0], "stale")

    def test_an_unreadable_record_is_said(self):
        with tempfile.TemporaryDirectory() as td:
            self.write_record(td, "{not json")
            self.assertEqual(sp.report_data(td)["md"], sp.REPORT_WORDS["unreadable"])
            self.assertEqual(sp.manifest_line(td, False)[0], "SKIPPED")

    def test_a_job_of_an_earlier_plan_never_writes_the_record(self):
        with tempfile.TemporaryDirectory() as td:
            plan = dict(plan_for(work=os.path.join(td, sp.WORK)), out=td, plan_id="old")
            with open(sp.record_path(td), "w") as fh:
                json.dump(dict(sp.skipped_record("turned off", "user")), fh)
            with mock.patch("sys.stderr", new=io.StringIO()):
                self.assertFalse(sp.write_owned(plan, dict(OK_REC)))
            self.assertEqual(sp.load_record(td)["status"], "skipped")
            with open(sp.record_path(td), "w") as fh:
                json.dump({"schema": sp.SCHEMA, "status": "planned", "plan_id": "old"}, fh)
            self.assertTrue(sp.write_owned(plan, dict(OK_REC)))
            self.assertEqual(sp.load_record(td)["status"], "ok")

    def test_skip_cancels_the_jobs_and_clears_the_last_run(self):
        with tempfile.TemporaryDirectory() as td:
            os.makedirs(os.path.join(td, sp.WORK))
            with open(os.path.join(td, sp.WORK, "jobs.txt"), "w") as fh:
                fh.write("build\t101\nimport\t102\njoin\t103\n")
            with open(sp.record_path(td), "w") as fh:
                json.dump({"schema": sp.SCHEMA, "status": "planned", "plan_id": "p"}, fh)
            for f in ("peaks_G.png", "T_top_hits.sky.zip", "notes.txt"):
                with open(os.path.join(td, f), "w") as fh:
                    fh.write("x")
            calls = []
            me = getpass.getuser()
            sacct = (f"101|{me}|sky_build|COMPLETED\n102_1|{me}|sky_import|RUNNING\n"
                     f"102_[2-3%3]|{me}|sky_import|PENDING\n103|{me}|sky_join|PENDING\n")

            def fake(argv, **k):
                calls.append(argv)
                return mock.Mock(returncode=0, stdout=sacct if argv[0] == "sacct" else "",
                                 stderr="")
            with mock.patch.object(sp.subprocess, "run", fake), \
                    mock.patch("sys.stdout", new=io.StringIO()), \
                    mock.patch("sys.stderr", new=io.StringIO()):
                sp.main(["skip", "--out", td, "--reason", "not wanted"])
            # only the live ones, after checking they are this user's and this step's
            self.assertEqual(calls[0][0], "sacct")
            self.assertEqual(calls[1], ["scancel", "102", "103"])
            self.assertFalse(os.path.exists(os.path.join(td, "peaks_G.png")))
            self.assertFalse(os.path.exists(os.path.join(td, "T_top_hits.sky.zip")))
            self.assertTrue(os.path.exists(os.path.join(td, "notes.txt")))   # not its own: kept
            self.assertFalse(os.path.isdir(os.path.join(td, sp.WORK)))

    def test_jobs_that_cannot_be_checked_or_are_not_ours_are_left_alone(self):
        with tempfile.TemporaryDirectory() as td:
            os.makedirs(os.path.join(td, sp.WORK))
            with open(os.path.join(td, sp.WORK, "jobs.txt"), "w") as fh:
                fh.write("build\t201\n")
            with open(sp.record_path(td), "w") as fh:
                json.dump({"schema": sp.SCHEMA, "status": "planned", "plan_id": "p"}, fh)
            for answer in (OSError("no sacct"),
                           mock.Mock(returncode=0, stdout="201|someone|other_job|RUNNING\n",
                                     stderr="")):
                calls = []

                def fake(argv, **k):
                    calls.append(argv)
                    if isinstance(answer, Exception):
                        raise answer
                    return answer
                with mock.patch.object(sp.subprocess, "run", fake), \
                        mock.patch("sys.stdout", new=io.StringIO()):
                    with self.assertRaises(SystemExit) as e:
                        sp.main(["skip", "--out", td, "--reason", "x"])
                self.assertIn("were not cancelled", str(e.exception))
                self.assertNotIn("scancel", [c[0] for c in calls])
                self.assertEqual(sp.load_record(td)["status"], "planned")

    def test_a_plan_that_cannot_run_leaves_no_old_pictures(self):
        with tempfile.TemporaryDirectory() as td:
            out = os.path.join(td, "output", "skyline")
            os.makedirs(out)
            with open(os.path.join(out, "peaks_OLD.png"), "w") as fh:
                fh.write("x")
            with open(sp.record_path(out), "w") as fh:
                json.dump(dict(OK_REC), fh)
            with mock.patch("sys.stdout", new=io.StringIO()):
                rc = sp.main(["plan", "--search", os.path.join(td, "nosearch"), "--out", out])
            self.assertEqual(rc, 3)
            self.assertFalse(os.path.exists(os.path.join(out, "peaks_OLD.png")))
            self.assertTrue(sp.report_data(os.path.dirname(out))["md"].startswith("**Not run.**"))

    def test_planned_and_submitted_are_told_apart(self):
        with tempfile.TemporaryDirectory() as td:
            self.write_record(td, {"schema": sp.SCHEMA, "status": "planned", "plan_id": "p"})
            self.assertEqual(sp.record_state(td)[0], "planned")
            self.assertEqual(sp.report_data(td)["md"], sp.REPORT_WORDS["planned"])
            os.makedirs(os.path.join(td, "skyline", sp.WORK))
            with open(os.path.join(td, "skyline", sp.WORK, "jobs.txt"), "w") as fh:
                fh.write("build\t1\n")
            self.assertEqual(sp.record_state(td)[0], "submitted")
            self.assertEqual(sp.report_data(td)["md"], sp.REPORT_WORDS["submitted"])
            self.assertNotIn("make the report again", sp.report_data(td)["md"])

    def test_a_crash_in_the_join_is_recorded(self):
        with tempfile.TemporaryDirectory() as td:
            w = os.path.join(td, sp.WORK)
            for d in ("status", "logs"):
                os.makedirs(os.path.join(w, d))
            p = dict(plan_for(raw=["/d/r1.raw"], work=w), out=td, plan_id="p1")
            with open(os.path.join(w, "plan.json"), "w") as fh:
                json.dump(p, fh)
            with open(sp.record_path(td), "w") as fh:
                json.dump({"schema": sp.SCHEMA, "status": "planned", "plan_id": "p1"}, fh)
            with open(os.path.join(w, "status", "build.json"), "w") as fh:
                json.dump({"ok": True, "plan_id": "p1"}, fh)
            with open(os.path.join(w, "status", "import_1.json"), "w") as fh:
                json.dump({"rc": 0, "plan_id": "p1"}, fh)
            for f in (os.path.join(td, "skyline_precursors.tsv"),
                      os.path.join(w, "skyline_transitions.tsv"),
                      os.path.join(w, "skyline_chromatograms.tsv"), p["share_zip"]):
                with open(f, "w") as fh:
                    fh.write("x")

            def boom(*a, **k):
                raise KeyError("selection")
            with mock.patch.object(sp, "run_skyline", lambda *a, **k: (0, "log")), \
                    mock.patch.object(sp, "finish", boom):
                self.assertEqual(sp.main(["join", "--work", w]), 1)
            rec = sp.load_record(td)
            self.assertEqual(rec["status"], "skipped")
            self.assertIn("the join job failed: KeyError", rec["reason"])

    def test_a_status_file_of_another_plan_is_not_taken(self):
        with tempfile.TemporaryDirectory() as td:
            w = os.path.join(td, sp.WORK)
            for d in ("status", "logs"):
                os.makedirs(os.path.join(w, d))
            p = dict(plan_for(raw=["/d/r1.raw"], work=w), out=td, plan_id="new")
            with open(os.path.join(w, "plan.json"), "w") as fh:
                json.dump(p, fh)
            with open(sp.record_path(td), "w") as fh:
                json.dump({"schema": sp.SCHEMA, "status": "planned", "plan_id": "new"}, fh)
            with open(os.path.join(w, "status", "build.json"), "w") as fh:
                json.dump({"ok": True, "plan_id": "old"}, fh)       # an earlier plan's job
            sp.main(["join", "--work", w])
            self.assertIn("build job ended without a result", sp.load_record(td)["reason"])

    def test_labels_carry_the_modifications(self):
        k = sp.precursor_key
        self.assertEqual(sp.key_label(k("PEPM[+16]K", 2)), "PEPM[+16]K/2+")
        self.assertNotEqual(sp.key_label(k("PEPM[+16]K", 2)), sp.key_label(k("PEPMK", 2)))
        self.assertEqual(sp.key_label(k("(UniMod:1)M(UniMod:35)ACK", 3)), "M[+58]ACK/3+")

    def test_the_searchs_own_mod_masses_match(self):
        # GlyGly is in no table of the skill: the search's --var-mod gives its mass
        flags = sp.command_flags(CMD + " --var-mod UniMod:121,114.042927,K")
        st = sp.search_settings(flags)
        self.assertEqual(st["mod_masses"], {121: 114.042927})
        sp.use_search_mod_masses(st["mod_masses"])
        self.assertEqual(sp.precursor_key("PEPK[+114]R", 2),
                         sp.precursor_key("PEPK(UniMod:121)R", 2))

    def test_submit_runs_once_per_plan(self):
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            w = os.path.join(td, sp.WORK)
            os.makedirs(w)
            with open(os.path.join(w, "jobs.txt"), "w") as fh:
                fh.write("build\t1\n")
            t = sp.submit_script(plan_for(work=w))
            r = subprocess.run(["bash", "-c", t], capture_output=True, text=True,
                               env=dict(os.environ, PATH="/usr/bin:/bin"))
            self.assertEqual(r.returncode, 2)
            self.assertIn("already submitted", r.stderr)
            with open(os.path.join(w, "jobs.txt")) as fh:
                self.assertEqual(fh.read(), "build\t1\n")             # never emptied

    def test_a_job_of_an_earlier_plan_does_nothing(self):
        with tempfile.TemporaryDirectory() as td:
            w = os.path.join(td, sp.WORK)
            os.makedirs(os.path.join(w, "status"))
            with open(os.path.join(w, "plan.json"), "w") as fh:
                json.dump(dict(plan_for(work=w), out=td, plan_id="new"), fh)
            with mock.patch("sys.stderr", new=io.StringIO()), \
                    self.assertRaises(SystemExit) as e:
                sp.main(["join", "--work", w, "--plan-id", "old"])
            self.assertEqual(e.exception.code, 0)
            self.assertFalse(os.path.exists(sp.record_path(td)))
            t = sp.job_script(dict(plan_for(), plan_id="p 1"), "join", ("high", "a", "q"))
            self.assertIn("--plan-id 'p 1'", t)

    def test_submit_records_each_job_as_it_is_made(self):
        import subprocess
        t = sp.submit_script(plan_for())
        self.assertEqual(t.count('>> "$J"'), 3)
        self.assertLess(t.index("build"), t.index("sky_2_import"))
        self.assertEqual(subprocess.run(["bash", "-n"], input=t, text=True).returncode, 0)

    def _cancel_case(self, td, sacct_out, scancel_rc=0, status="planned"):
        os.makedirs(os.path.join(td, sp.WORK), exist_ok=True)
        with open(os.path.join(td, sp.WORK, "jobs.txt"), "w") as fh:
            fh.write("build\t301\nimport\t302\n")
        with open(sp.record_path(td), "w") as fh:
            json.dump({"schema": sp.SCHEMA, "status": status, "plan_id": "p",
                       "reason": "x"}, fh)
        calls = []

        def fake(argv, **k):
            calls.append(argv)
            if argv[0] == "sacct":
                return mock.Mock(returncode=0, stdout=sacct_out, stderr="")
            return mock.Mock(returncode=scancel_rc, stdout="", stderr="scancel: error")
        return calls, fake

    def test_a_failing_scancel_refuses(self):
        me = getpass.getuser()
        with tempfile.TemporaryDirectory() as td:
            calls, fake = self._cancel_case(td, f"301|{me}|sky_build|RUNNING\n302|{me}|sky_import|PENDING\n", 1)
            with mock.patch.object(sp.subprocess, "run", fake), \
                    mock.patch("sys.stdout", new=io.StringIO()), self.assertRaises(SystemExit) as e:
                sp.main(["skip", "--out", td, "--reason", "x"])
            self.assertIn("scancel 301 302 failed", str(e.exception))
            self.assertEqual(sp.load_record(td)["status"], "planned")

    def test_a_job_accounting_does_not_know_yet(self):
        me = getpass.getuser()
        with tempfile.TemporaryDirectory() as td:   # just submitted: refuse
            calls, fake = self._cancel_case(td, f"301|{me}|sky_build|RUNNING\n")
            with mock.patch.object(sp.subprocess, "run", fake), \
                    mock.patch("sys.stdout", new=io.StringIO()), self.assertRaises(SystemExit) as e:
                sp.main(["skip", "--out", td, "--reason", "x"])
            self.assertIn("does not know job(s) 302", str(e.exception))
        with tempfile.TemporaryDirectory() as td:   # a finished run long ago: gone is finished
            calls, fake = self._cancel_case(td, "", status="skipped")
            with mock.patch.object(sp.subprocess, "run", fake), \
                    mock.patch("sys.stdout", new=io.StringIO()), \
                    mock.patch("sys.stderr", new=io.StringIO()):
                sp.main(["skip", "--out", td, "--reason", "x"])
            self.assertNotIn("scancel", [c[0] for c in calls])
            self.assertIn("Turned off", sp.load_record(td)["reason"])

    def test_replan_while_jobs_run_cancels_them(self):
        me = getpass.getuser()
        with tempfile.TemporaryDirectory() as td:
            out = os.path.join(td, "output", "skyline")
            calls, fake = self._cancel_case(
                out, f"301|{me}|sky_build|COMPLETED\n302_[1-2]|{me}|sky_import|PENDING\n")
            with mock.patch.object(sp.subprocess, "run", fake), \
                    mock.patch("sys.stdout", new=io.StringIO()), \
                    mock.patch("sys.stderr", new=io.StringIO()):
                rc = sp.main(["plan", "--search", os.path.join(td, "nosearch"), "--out", out,
                              "--replan"])
            self.assertEqual(rc, 3)                       # this plan cannot run: no search
            self.assertIn(["scancel", "302"], calls)
            self.assertFalse(os.path.isdir(os.path.join(out, sp.WORK)))

    def test_the_record_says_not_finished_before_the_old_run_is_removed(self):
        with tempfile.TemporaryDirectory() as td:
            out = os.path.join(td, "output", "skyline")
            os.makedirs(out)
            with open(sp.record_path(out), "w") as fh:
                json.dump(dict(OK_REC), fh)
            seen = []

            def clear(o):
                seen.append(sp.load_record(o)["status"])
                raise PermissionError("cannot remove")
            with mock.patch.object(sp, "clear_previous", clear), self.assertRaises(PermissionError):
                sp.main(["plan", "--search", os.path.join(td, "s"), "--out", out])
            self.assertEqual(seen, ["skipped"])               # not the old "ok"
            self.assertIn("PermissionError", sp.load_record(out)["reason"])

    def test_a_plan_stopped_by_an_error_says_so(self):
        with tempfile.TemporaryDirectory() as td:
            search = TestNotRun().make_search(td)
            out = os.path.join(td, "output", "skyline")
            with mock.patch.object(sp, "make_plan", side_effect=SystemExit("bad --account")), \
                    self.assertRaises(SystemExit):
                sp.main(["plan", "--search", search, "--out", out])
            rec = sp.load_record(out)
            self.assertEqual(rec["status"], "skipped")
            self.assertIn("bad --account", rec["reason"])

    def test_a_failed_join_removes_its_half_made_outputs(self):
        with tempfile.TemporaryDirectory() as td:
            w = os.path.join(td, sp.WORK)
            for d in ("status", "logs"):
                os.makedirs(os.path.join(w, d))
            p = dict(plan_for(raw=["/d/r1.raw"], work=w), out=td, plan_id="p1")
            with open(os.path.join(w, "plan.json"), "w") as fh:
                json.dump(p, fh)
            with open(sp.record_path(td), "w") as fh:
                json.dump({"schema": sp.SCHEMA, "status": "planned", "plan_id": "p1"}, fh)
            for name, body in (("build.json", {"ok": True, "plan_id": "p1"}),
                               ("import_1.json", {"rc": 0, "plan_id": "p1"})):
                with open(os.path.join(w, "status", name), "w") as fh:
                    json.dump(body, fh)
            for f in (os.path.join(td, "skyline_precursors.tsv"), os.path.join(td, "peaks_G.png"),
                      os.path.join(w, "skyline_transitions.tsv"),
                      os.path.join(w, "skyline_chromatograms.tsv"), p["share_zip"]):
                with open(f, "w") as fh:
                    fh.write("x")
            with mock.patch.object(sp, "run_skyline", lambda *a, **k: (0, "log")), \
                    mock.patch.object(sp, "finish", side_effect=ValueError("bad table")):
                self.assertEqual(sp.main(["join", "--work", w, "--plan-id", "p1"]), 1)
            left = [f for f in os.listdir(td) if sp.OWN_OUTPUTS.match(f)]
            self.assertEqual(left, [])
            self.assertIn("ValueError: bad table", sp.load_record(td)["reason"])

    def test_the_disagreement_picture_shows_the_runs_that_disagreed(self):
        runs = [f"run{i:02d}" for i in range(30)]
        flagged = [{"precursor": "X/2+", "run": "run29"}, {"precursor": "X/2+", "run": "run27"},
                   {"precursor": "Y/3+", "run": "run28"}]
        got = sp.disagreement_runs(flagged, {"X/2+"}, runs, {})
        self.assertEqual(len(got), sp.MAX_RUNS_DRAWN)
        self.assertIn("run29", got)
        self.assertIn("run27", got)
        self.assertNotIn("run28", got)               # another precursor's disagreement

    def test_dda_from_the_manifest(self):
        with tempfile.TemporaryDirectory() as td:
            search = TestNotRun().make_search(td)
            man = os.path.join(td, "workflow.manifest.json")
            with open(man, "w") as fh:
                json.dump({"acquisition": "DDA", "instrument_class": "orbitrap_generic"}, fh)
            out = os.path.join(td, "output", "skyline")
            with mock.patch("sys.stdout", new=io.StringIO()):
                self.assertEqual(sp.main(["plan", "--search", search, "--out", out,
                                          "--manifest", man]), 3)
            self.assertIn("DDA", sp.load_record(out)["reason"])

    def test_nothing_matched_is_not_ok(self):
        with tempfile.TemporaryDirectory() as td:
            w = os.path.join(td, sp.WORK)
            os.makedirs(w)
            with open(os.path.join(td, "skyline_precursors.tsv"), "w", newline="") as fh:
                wr = csv.DictWriter(fh, SKY_HEADER, delimiter="\t")
                wr.writeheader()
                wr.writerow(sky_row("sp|P1|A_HUMAN", "PEPTIDEK", 2, "other_run", 10, 9.9, 10.1))
            with open(os.path.join(w, "diann_targets.tsv"), "w", newline="") as fh:
                wr = csv.DictWriter(fh, list(dia_row("P1", "PEPTIDEK", 2, "r1", 1, 1, 1)),
                                    delimiter="\t")
                wr.writeheader()
                wr.writerow(dia_row("P1", "PEPTIDEK", 2, "r1", 10.0, 9.9, 10.1))
            p = dict(plan_for(work=w), out=td)
            build = {"selection": {"targets": [{"protein_group": "P1", "gene": "G1"}]}}
            with self.assertRaises(sp.Skip) as e:
                sp.finish(p, build, ["/d/r1.raw"], [])
            self.assertIn("nothing could be compared", str(e.exception))
            with open(os.path.join(w, "diann_targets.tsv"), "w") as fh:   # DIA-NN: no rows
                fh.write("Run\tProtein.Group\tModified.Sequence\tPrecursor.Charge\n")
            with self.assertRaises(sp.Skip):
                sp.finish(p, build, ["/d/r1.raw"], [])

    def test_what_could_not_be_compared_is_said(self):
        md = sp.report_section(dict(OK_REC, results=dict(OK_REC["results"], not_compared=7)))["md"]
        self.assertIn("Another 7 of Skyline's precursor-runs had no DIA-NN identification", md)

    def test_the_median_is_over_the_rows_skyline_agreed_on(self):
        # the reviewer's numbers: 1, 2, 50, 60 and 70 s -> 50, 60, 70 flagged; median of 1, 2
        rows = [{"apex_difference_s": d, "reason": "" if d <= sp.APEX_DIFF_S else "apex far"}
                for d in (1, 2, 50, 60, 70)]
        s = sp.agreement_summary(rows)
        self.assertEqual(s["apex_median_abs_s"], 1.5)
        self.assertEqual(s["within_s_pct"], 40.0)
        md = sp.report_section(dict(OK_REC))["md"]
        self.assertIn("where it picked the same peak, the median difference", md)


OPEN_WORDS = ("Skyline 26.1 or later on Windows", "https://skyline.ms",
              "double-click it, or use File › Open, and Skyline unpacks it",
              "the raw files are not needed to look at the peaks",
              "needed only to re-pick peaks or add peptides")


class TestOpeningTheDocument(unittest.TestCase):
    """Brett: the Skyline file easy to find, and openable by people with Skyline on Windows.
    SkylineCmd 26.1's --help has no option to write an older document format, so the oldest
    desktop Skyline is 26.1."""

    def test_the_version_is_the_pinned_release(self):
        self.assertEqual(sp.desktop_version(), "26.1")

    def test_the_report_names_the_file_by_its_delivered_path(self):
        md = sp.report_section(dict(OK_REC))["md"]
        self.assertIn("The Skyline document is `skyline/T_top_hits.sky.zip`.", md)
        for w in OPEN_WORDS:
            self.assertIn(w, md)
        self.assertEqual(client_text.peak_check_overclaims(md), [])

    def test_output_files_md(self):
        import make_report
        desc = make_report.describe("T_top_hits.sky.zip")[1]
        for w in OPEN_WORDS:
            self.assertIn(w, desc)

    def test_the_delivery_readme(self):
        import core_submission as cs
        from test_core_submission import write_summary
        with tempfile.TemporaryDirectory() as td:
            _p, summary = write_summary(td)
            delivered = {"Analysis_Report.html", "figures/pca.png", "skyline/peaks_G.png",
                         "skyline/T_top_hits.sky.zip"}
            for for_html in (False, True):
                text = cs.build_readme(summary, delivered, [], "analysis", for_html=for_html)
                self.assertIn("skyline/T_top_hits.sky.zip", text)
                for w in OPEN_WORDS:
                    self.assertIn(w, text)
            self.assertNotIn("Skyline", cs.build_readme(summary, {"Analysis_Report.html"}, [],
                                                         "analysis"))

    def test_the_session_readme(self):
        import session_docs
        from test_session_docs import Agents, DPC
        with tempfile.TemporaryDirectory() as td:
            p = Agents().build(td, DPC)
            sky = os.path.join(p["output_dir"], "skyline")
            os.makedirs(sky)
            open(os.path.join(sky, "T_top_hits.sky.zip"), "w").close()
            sp.write_json(sp.record_path(sky), dict(OK_REC))       # a finished check (readme_rows)
            md = session_docs.readme_md(session_docs.gather(p["session_dir"]))
            self.assertIn("| `output/skyline/T_top_hits.sky.zip` | the Skyline document", md)
            for w in OPEN_WORDS:
                self.assertIn(w, md)
            # a check that did not finish: one row with its state, never the pictures
            sp.write_json(sp.record_path(sky), sp.skipped_record("the container is gone", "build"))
            md = session_docs.readme_md(session_docs.gather(p["session_dir"]))
        self.assertIn("| `output/skyline` | the peak check (Skyline): not run: the container is "
                      "gone |", md)
        self.assertNotIn("T_top_hits.sky.zip", md)


class TestReadmeRows(unittest.TestCase):
    """2.14.2: readme_rows() described any output/skyline/ folder as "pictures of the checked
    proteins' peaks ... and the Skyline document", even one holding only a skip record. The rows
    now follow record_state(), the rule the report, the Methods and the delivery read: only a
    finished check is described as pictures and a document; any other state gets ONE row that
    says what state the check is in and why (rule 1: the record describes itself)."""

    def rows(self, rec, files=()):
        with tempfile.TemporaryDirectory() as td:
            sky = os.path.join(td, "skyline")
            os.makedirs(sky)
            for f in files:
                open(os.path.join(sky, f), "w").close()
            if rec is not None:
                sp.write_json(sp.record_path(sky), rec)
            return sp.readme_rows(sky, "output/skyline")

    def test_finished(self):
        self.assertEqual(self.rows(dict(OK_REC), ["T_top_hits.sky.zip", "peaks_G1.png"]), [
            ("output/skyline", "the peak check: pictures of the checked proteins' peaks, drawn "
                               "with Skyline, its tables and the Skyline document"),
            ("output/skyline/T_top_hits.sky.zip", "the Skyline document of the checked proteins. "
             + sp.open_document_text())])
        # finished, but no document in the folder: the row does not promise one
        self.assertEqual(self.rows(dict(OK_REC), ["peaks_G1.png"]), [
            ("output/skyline", "the peak check: pictures of the checked proteins' peaks, drawn "
                               "with Skyline, and its tables")])

    def test_skipped(self):
        rows = self.rows(sp.skipped_record("Turned off for this analysis: tables only", "user"),
                         ["T_top_hits.sky.zip"])
        self.assertEqual(rows, [("output/skyline", "the peak check (Skyline): not run: Turned "
                                                   "off for this analysis: tables only")])

    def test_failed(self):
        # a job that failed writes a skipped record with the stage and the reason (_fail)
        why = "Skyline could not build the document (exit 2: Error: no library)"
        rows = self.rows(sp.skipped_record(why, "build"), ["T_top_hits.sky.zip", "peaks_G.png"])
        self.assertEqual(rows, [("output/skyline", f"the peak check (Skyline): not run: {why}")])
        # a record whose status is no state of the check's ("failed") cannot be read as one
        rows = self.rows(dict(OK_REC, status="failed"))
        self.assertEqual(rows, [("output/skyline", "the peak check (Skyline): "
                                 + sp.STATE_WORDS["unreadable"])])

    def test_the_other_states(self):
        planned = {"schema": sp.SCHEMA, "status": "planned", "plan_id": "p"}
        self.assertEqual(self.rows(planned), [("output/skyline", "the peak check (Skyline): not "
                                               "finished: " + sp.STATE_WORDS["planned"])])
        rows = self.rows(None, ["peaks_G.png"])
        self.assertEqual(len(rows), 1)
        self.assertIn(sp.RECORD, rows[0][1])
        self.assertNotIn("pictures of the checked", rows[0][1])
        with tempfile.TemporaryDirectory() as td:
            self.assertEqual(sp.readme_rows(os.path.join(td, "skyline"), "output/skyline"), [])

    def test_the_rows_follow_the_delivery_rule(self):
        """L9 review F4: readme_rows decides with deliverable() -- the rule the delivery ships
        skyline/ by -- not with a record_state() test of its own, so a change to the rule can
        never leave the README describing a folder the delivery does not ship."""
        with tempfile.TemporaryDirectory() as td:
            sky = os.path.join(td, "skyline")
            os.makedirs(sky)
            open(os.path.join(sky, "T_top_hits.sky.zip"), "w").close()
            sp.write_json(sp.record_path(sky), dict(OK_REC))
            with mock.patch.object(sp, "deliverable",
                                   return_value=(False, "stale", "not shipped")) as d:
                rows = sp.readme_rows(sky, "output/skyline")
            d.assert_called_once_with(os.path.abspath(td))
        self.assertEqual(len(rows), 1)
        self.assertNotIn("pictures of the checked", rows[0][1])

    def test_one_branch_from_a_state_to_its_words(self):
        """L9 review F4 (rule 3): the README row and deliverable()'s `why` are state_reason()'s
        two framings of ONE branch; deliverable() has no branch of its own."""
        S = sp.SCHEMA
        cases = {"missing": None, "unreadable": {"schema": S, "status": "weird"},
                 "skipped": {"schema": S, "status": "skipped", "reason": "turned off"},
                 "planned": {"schema": S, "status": "planned", "plan_id": "p"},
                 "stale": dict(OK_REC, de_record_sha256="0" * 64)}
        for state, rec in cases.items():
            with self.subTest(state=state), tempfile.TemporaryDirectory() as td:
                sky = os.path.join(td, "skyline")
                os.makedirs(sky)
                if rec is not None:
                    sp.write_json(sp.record_path(sky), rec)
                got_state, got_rec = sp.record_state(td)
                self.assertEqual(got_state, state)
                self.assertEqual(sp.deliverable(td),
                                 (False, state, sp.state_reason(state, got_rec, standalone=True)))
                self.assertEqual(sp.readme_rows(sky, "output/skyline"),
                                 [("output/skyline", "the peak check (Skyline): "
                                   + sp.state_reason(state, got_rec))])
        with open(os.path.splitext(sp.__file__)[0] + ".py", encoding="utf-8") as fh:
            src = fh.read()
        self.assertEqual(src.count("STATE_WORDS.get("), 1, "one branch for the words")


class TestOneFolderSentence(unittest.TestCase):
    """L9 review F5 (rule 3): the README sentence describing a finished check's folder, and its
    Skyline document's row, are defined once (skyline_peaks.folder_rows); the delivery README
    reads them instead of carrying its own copy."""

    def test_the_delivery_readme_reads_the_one_definition(self):
        import core_submission as cs
        from test_core_submission import write_summary
        with tempfile.TemporaryDirectory() as td:
            _p, summary = write_summary(td)
            delivered = {"Analysis_Report.html", "figures/pca.png", "skyline/peaks_G.png",
                         "skyline/T_top_hits.sky.zip"}
            with mock.patch.object(sp, "folder_rows",
                                   return_value=[("skyline/", "ONE DEFINITION")]) as fr:
                text = cs.build_readme(summary, delivered, [], "analysis")
            fr.assert_called_once_with("skyline/", ["skyline/T_top_hits.sky.zip"])
            self.assertIn("ONE DEFINITION", text)
            self.assertNotIn("pictures of the checked", text)
            # and the real rows, word for word, with or without a document
            for docs in (["skyline/T_top_hits.sky.zip"], []):
                names = {"Analysis_Report.html", "figures/pca.png", "skyline/peaks_G.png"} | set(docs)
                text = cs.build_readme(summary, names, [], "analysis")
                for _path, words in sp.folder_rows("skyline/", docs):
                    self.assertIn(f" | {words} |", text)
        with open(os.path.splitext(cs.__file__)[0] + ".py", encoding="utf-8") as fh:
            self.assertNotIn("pictures of the checked", fh.read())

    def test_the_session_readme_reads_it_too(self):
        with tempfile.TemporaryDirectory() as td:
            sky = os.path.join(td, "skyline")
            os.makedirs(sky)
            open(os.path.join(sky, "T_top_hits.sky.zip"), "w").close()
            sp.write_json(sp.record_path(sky), dict(OK_REC))
            self.assertEqual(sp.readme_rows(sky, "output/skyline"),
                             sp.folder_rows("output/skyline",
                                            ["output/skyline/T_top_hits.sky.zip"]))


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


class TestMethods(unittest.TestCase):
    def test_the_paragraph_only_when_it_ran(self):
        import subprocess
        for status in ("ok", "skipped"):
            with tempfile.TemporaryDirectory() as td:
                tables, sky = os.path.join(td, "output", "tables"), os.path.join(td, "output",
                                                                                "skyline")
                os.makedirs(tables)
                os.makedirs(sky)
                with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
                    json.dump({"adjp": 0.05, "method": "dpc"}, fh)
                rec = dict(OK_REC, status=status, reason="turned off")
                rec["methods_paragraph"] = sp.methods_paragraph(dict(OK_REC))
                with open(os.path.join(sky, sp.RECORD), "w") as fh:
                    json.dump(rec, fh)
                raw = os.path.join(td, "FL_sample01.raw")
                open(raw, "w").close()
                out = os.path.join(td, "methods.md")
                r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                                    "--raw", raw, "--de-dir", tables, "--out", out],
                                   capture_output=True, text=True, timeout=300)
                self.assertEqual(r.returncode, 0, r.stderr)
                with open(out, encoding="utf-8") as fh:
                    text = fh.read()
                if status == "ok":
                    self.assertIn("## Peak check (Skyline)", text)
                    self.assertIn("not an independent identification", text)
                    self.assertIn(SKYLINE_PIN["docker_digest"], text)
                else:
                    self.assertNotIn("Peak check", text)


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

    def test_a_second_request_replaces_the_first(self):
        with tempfile.TemporaryDirectory() as td:
            w = os.path.join(td, "skyline", sp.WORK)
            dest = os.path.join(td, sp.FOLLOWUP_DIR)
            for d in (os.path.join(w, "prm"), os.path.join(w, "logs"), dest):
                os.makedirs(d)
            p = dict(plan_for(work=w), out=os.path.join(td, "skyline"))
            with open(os.path.join(w, "plan.json"), "w") as fh:
                json.dump(p, fh)
            targets = sp.prm_plan({"family": "orbitrap"}, ["thermo_fusion"], True, None)
            with open(os.path.join(w, "prm", "prm.json"), "w") as fh:
                json.dump({"targets": targets, "dest": dest}, fh)
            for f in (os.path.join(w, "prm", "refined.sky"), os.path.join(dest, targets[0]["file"])):
                with open(f, "w") as fh:
                    fh.write("old")

            def fake(plan, args, log_name, rw=(), timeout_s=None):
                # Skyline refuses an --out that exists: the last request's files are gone first
                self.assertFalse(os.path.exists(os.path.join(w, "prm", "refined.sky")))
                self.assertEqual(rw, [dest])
                with open(os.path.join(dest, targets[0]["file"]), "w") as fh:
                    fh.write("Compound,m/z\nA,1\n")
                return 0, os.path.join(w, "logs", log_name)
            with mock.patch.object(sp, "run_skyline", fake):
                self.assertEqual(sp.main(["prm-run", "--work", w]), 0)
            with open(os.path.join(dest, "followup_assay.json")) as fh:
                self.assertEqual(json.load(fh)["rows"], {targets[0]["file"]: 1})

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


class TestImportTimeAndDocument(unittest.TestCase):
    """2.14.2, I39. On HIVE's high, import tasks suspended by Slurm for ~28 min were killed by the
    skill's own wall-clock timer (rc 124) seconds after they resumed, because suspended time
    counted against it: Slurm's --time (which does not count it) now ends a job, and only
    `SkylineCmd --version` has a time limit of its own. And a task run again after the join
    (by hand, or requeued) copied the document the join had saved its results into: each task
    now copies the build's pristine copy, <stem>.built.sky."""

    def work(self, td):
        w = os.path.join(td, sp.WORK)
        for d in ("status", "logs"):
            os.makedirs(os.path.join(w, d))
        p = dict(plan_for(work=w), out=td, plan_id="p1", fasta="/db/search.fasta",
                 report="/s/report.parquet", q_cutoff=0.01)
        with open(os.path.join(w, "plan.json"), "w") as fh:
            json.dump(p, fh)
        return w, p

    @staticmethod
    def put(path, text):
        with open(path, "w") as fh:
            fh.write(text)

    def test_import_call_has_no_timeout(self):
        with tempfile.TemporaryDirectory() as td:
            w, p = self.work(td)
            stem = os.path.splitext(p["doc"])[0]
            self.put(stem + ".sky", "<srm_settings/>")
            self.put(stem + sp.BUILT_TAG + ".sky", "<srm_settings/>")
            sif = os.path.join(td, "other.sif")             # not the Core's: no sha256 to match
            self.put(sif, "image")
            calls = []

            class FakePopen:               # SkylineCmd, started in a process group of its own
                def __init__(self, argv, stdout=None, stderr=None, **k):
                    self.argv, self.k, self.pid = argv, k, 0
                    if "--version" in argv:
                        stdout.write("Skyline (64-bit) 26.1.0.266 (578b175)\n")
                    for a in argv:
                        if a.startswith("--in="):                 # the task's own copy
                            tdir = os.path.dirname(a.split("=", 1)[1][2:].replace("\\", "/"))
                            open(os.path.join(tdir, "r1_1.raw.skyd"), "w").close()

                def wait(self, timeout=None):
                    calls.append((self.argv, dict(self.k, timeout=timeout)))
                    return 0
            with mock.patch.object(sp.subprocess, "Popen", FakePopen), \
                    mock.patch("sys.stderr", new=io.StringIO()):
                self.assertEqual(sp.main(["import", "--work", w, "--task", "1",
                                          "--plan-id", "p1"]), 0)
                sp.run_skyline(p, ["--in=" + sp.win(p["doc"])], "build.txt")
                sp.verify_container(dict(p, sif=sif))
            imp, other, version = calls
            self.assertTrue(any(a.startswith("--import-file=") for a in imp[0]))
            self.assertIsNone(imp[1].get("timeout"))
            self.assertIsNone(other[1].get("timeout"))
            self.assertIn("--version", version[0])
            self.assertEqual(version[1].get("timeout"), sp.VERSION_TIMEOUT_S)
            self.assertEqual(sp.VERSION_TIMEOUT_S, 600)
            for c in calls:                # its own group: end_child ends all of it (L9 F1)
                self.assertTrue(c[1].get("start_new_session"))

    def test_import_copies_the_pristine_document_not_the_joined_one(self):
        with tempfile.TemporaryDirectory() as td:
            w, p = self.work(td)
            stem = os.path.splitext(p["doc"])[0]
            for ext, joined, built in ((".sky", "<joined/>", "<built/>"),
                                       (".blib", "joined lib", "built lib")):
                self.put(stem + ext, joined)
                self.put(stem + sp.BUILT_TAG + ext, built)
            seen = {}

            def fake(plan, args, log_name, rw=(), timeout_s=None):
                doc = args[0].split("=", 1)[1][2:].replace("\\", "/")
                tdir = os.path.dirname(doc)
                seen.update({f: open(os.path.join(tdir, f)).read() for f in os.listdir(tdir)})
                open(os.path.join(tdir, "r2_1.raw.skyd"), "w").close()
                return 0, os.path.join(w, "logs", log_name)
            with mock.patch.object(sp, "run_skyline", fake):
                self.assertEqual(sp.main(["import", "--work", w, "--task", "2"]), 0)
            # under the document's own names: Skyline finds the library by the document's name
            self.assertEqual(seen, {"T_top_hits.sky": "<built/>", "T_top_hits.blib": "built lib"})
            self.assertEqual(sp.pristine_doc(p), stem + sp.BUILT_TAG + ".sky")
            with open(os.path.join(w, "status", "import_2.json")) as fh:
                self.assertEqual(json.load(fh)["document"], stem + sp.BUILT_TAG + ".sky")

    def test_old_plan_without_pristine_copy_warns(self):
        with tempfile.TemporaryDirectory() as td:
            w, p = self.work(td)
            stem = os.path.splitext(p["doc"])[0]
            self.put(stem + ".sky", "<only/>")
            seen = []

            def fake(plan, args, log_name, rw=(), timeout_s=None):
                doc = args[0].split("=", 1)[1][2:].replace("\\", "/")
                seen.append(open(doc).read())
                open(os.path.join(os.path.dirname(doc), "r1_1.raw.skyd"), "w").close()
                return 0, os.path.join(w, "logs", log_name)
            with mock.patch.object(sp, "run_skyline", fake), \
                    mock.patch("sys.stderr", new=io.StringIO()) as err:
                self.assertEqual(sp.main(["import", "--work", w, "--task", "1"]), 0)
            self.assertEqual(seen, ["<only/>"])
            self.assertIn("WARNING: no pristine copy; importing from the document the join may "
                          "have saved over", err.getvalue())
            with open(os.path.join(w, "status", "import_1.json")) as fh:
                self.assertEqual(json.load(fh)["document"], p["doc"])

    def test_the_build_keeps_a_pristine_copy(self):
        with tempfile.TemporaryDirectory() as td:
            w, p = self.work(td)
            stem = os.path.splitext(p["doc"])[0]

            def fake(plan, args, log_name, rw=(), timeout_s=None):
                self.assertIsNone(timeout_s)
                self.put(stem + ".sky", '<srm_settings><peptide sequence="PEPTIDEK"/>')
                self.put(stem + ".blib", "lib")
                return 0, os.path.join(w, "logs", log_name)
            sel = {"mode": "de_hits", "targets": [{"protein_group": "P1", "gene": "G1"}]}
            with mock.patch.object(sp, "needs_for_pictures", lambda: None), \
                    mock.patch.object(sp, "verify_container", lambda plan: {"name": "Skyline"}), \
                    mock.patch.object(sp, "select_targets", lambda plan: dict(sel)), \
                    mock.patch.object(sp, "fasta_subset", lambda *a: {"P1"}), \
                    mock.patch.object(sp, "diann_rows", lambda *a: (4, 30.0)), \
                    mock.patch.object(sp, "run_skyline", fake), \
                    mock.patch("sys.stderr", new=io.StringIO()):
                self.assertEqual(sp.main(["build", "--work", w, "--plan-id", "p1"]), 0)
            for ext in (".sky", ".blib"):
                with open(stem + sp.BUILT_TAG + ext) as a, open(stem + ext) as b:
                    self.assertEqual(a.read(), b.read())
            with open(os.path.join(w, "status", "build.json")) as fh:
                self.assertTrue(json.load(fh)["ok"])

    def test_a_join_run_twice_starts_from_the_pristine_document(self):
        """L9 review F6: the join saves the replicates into the document (--save). Run again --
        requeued on low after its save, or by hand -- it must not re-open that joined document
        and import them a second time: it restores the build's pristine copy first, and the
        joined cache goes with it; the per-file caches stay."""
        with tempfile.TemporaryDirectory() as td:
            w, p = self.work(td)
            p["raw_files"] = ["/d/r1.raw"]
            with open(os.path.join(w, "plan.json"), "w") as fh:
                json.dump(p, fh)
            stem = os.path.splitext(p["doc"])[0]
            for ext, text in ((".sky", "<built/>"), (".blib", "built lib")):
                self.put(stem + ext, text)
                self.put(stem + sp.BUILT_TAG + ext, text)
            self.put(os.path.join(w, "r1.raw.skyd"), "the import's per-file cache")
            sp.write_json(sp.record_path(td), {"schema": sp.SCHEMA, "status": "planned",
                                               "plan_id": "p1"})
            for name, body in (("build.json", {"ok": True, "plan_id": "p1"}),
                               ("import_1.json", {"rc": 0, "plan_id": "p1"})):
                sp.write_json(os.path.join(w, "status", name), body)
            seen = []

            def fake(plan, args, log_name, rw=(), timeout_s=None):
                # what the batch's --in opens, as the join's SkylineCmd call finds it
                seen.append((open(p["doc"]).read(), os.path.exists(stem + ".skyd"),
                             os.path.exists(stem + ".sky.view"),
                             os.path.exists(os.path.join(w, "r1.raw.skyd"))))
                self.put(p["doc"], "<joined/>")                  # its --save
                self.put(stem + ".skyd", "joined cache")
                self.put(stem + ".sky.view", "view")
                for f in (os.path.join(td, "skyline_precursors.tsv"),
                          os.path.join(w, "skyline_transitions.tsv"),
                          os.path.join(w, "skyline_chromatograms.tsv"), p["share_zip"]):
                    self.put(f, "x")
                return 0, os.path.join(w, "logs", log_name)
            done = {"schema": sp.SCHEMA, "status": "ok", "results": {"flagged": 0, "compared": 1}}
            with mock.patch.object(sp, "run_skyline", fake), \
                    mock.patch.object(sp, "finish", lambda *a: dict(done)), \
                    mock.patch("sys.stderr", new=io.StringIO()):
                for _ in range(2):
                    self.assertEqual(sp.main(["join", "--work", w, "--plan-id", "p1"]), 0)
            self.assertEqual(seen, [("<built/>", False, False, True)] * 2)
            self.assertEqual(open(stem + sp.BUILT_TAG + ".sky").read(), "<built/>")
            self.assertEqual(sp.load_record(td)["status"], "ok")


# A stand-in for `apptainer exec ... wine SkylineCmd <args>` (TestStoppedBySlurm), found on PATH
# by the real job: it writes what a half-finished SkylineCmd would have left (FAKE_SKY_WRITES),
# starts a grandchild that keeps writing (FAKE_SKY_BEAT; Wine's processes), says both pids
# (FAKE_SKY_PIDS), sends the job SIGTERM itself when FAKE_SKY_KILL_PARENT is set, and hangs.
FAKE_APPTAINER = """#!{python}
import os, signal, subprocess, sys, time
if "--version" in sys.argv:
    print("Skyline (64-bit) 26.1.0.266 (578b175)")
    sys.exit(0)
for p in filter(None, os.environ.get("FAKE_SKY_WRITES", "").split(os.pathsep)):
    with open(p, "w") as fh:
        fh.write("half")
beat = subprocess.Popen([sys.executable, "-c", "import sys, time\\n"
                         "while True:\\n    open(sys.argv[1], 'a').write('.')\\n    time.sleep(0.1)",
                         os.environ["FAKE_SKY_BEAT"]])
while not os.path.exists(os.environ["FAKE_SKY_BEAT"]):     # it writes before anyone looks
    time.sleep(0.01)
pids = os.environ["FAKE_SKY_PIDS"]
with open(pids + ".tmp", "w") as fh:
    fh.write(f"{os.getpid()} {beat.pid}")
os.replace(pids + ".tmp", pids)
if os.environ.get("FAKE_SKY_KILL_PARENT"):
    os.kill(os.getppid(), signal.SIGTERM)
time.sleep(300)
"""


# A stand-in for `scontrol show job -o <id>` inside a job: HIVE's form of the line
# ("... RunTime=00:00:05 TimeLimit=00:15:00 TimeMin=N/A ...", h2 job 25387830), with the limit
# the test gives (FAKE_SCONTROL_LIMIT).
FAKE_SCONTROL = """#!{python}
import os, sys
assert sys.argv[1:4] == ["show", "job", "-o"], sys.argv
print(f"JobId={sys.argv[4]} JobName=sky RunTime=00:00:05 "
      f"TimeLimit={os.environ['FAKE_SCONTROL_LIMIT']} TimeMin=N/A")
"""


class TestStoppedBySlurm(unittest.TestCase):
    """L9 review F1, F3, F7. With no timer of the skill's own (I39), Slurm's --time is what ends a
    hung job: SIGTERM, then SIGKILL. Python's default SIGTERM action ended the job on the spot --
    no except path, no record -- so a join stopped mid-batch left the record "submitted" for good
    (the report: "make the report again once they finish") beside whatever SkylineCmd had
    already written, and finalize zipped it. Each job now handles SIGTERM and the early warning
    (#SBATCH --signal=B:USR1@120) like any failure: SkylineCmd's process group is ended, the
    partial outputs go, and the record (an import: its status file; a follow-up: its README)
    says the job was stopped at its time limit.

    The job runs as the real batch job does, as its own process (the job script exec's it), with
    a stand-in apptainer on PATH that writes partial outputs, starts a grandchild and hangs."""

    def setUp(self):
        self.pids = []
        self.td = tempfile.mkdtemp()
        self.addCleanup(__import__("shutil").rmtree, self.td, True)

    def tearDown(self):
        import signal
        for pid in self.pids:                    # nothing of a failed run is left running
            try:
                os.kill(pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass

    def run_job(self, argv, writes, kill_parent=False, send=None, limit=None):
        """The job as its own process. `limit`: the TimeLimit a stand-in scontrol reports for
        it (and SLURM_JOB_ID set); None: not in a job, so no limit can be read."""
        import signal
        import subprocess
        import time
        b = os.path.join(self.td, "bin")
        os.makedirs(b, exist_ok=True)
        fakes = {"apptainer": FAKE_APPTAINER}
        if limit:
            fakes["scontrol"] = FAKE_SCONTROL
        for name, text in fakes.items():
            with open(os.path.join(b, name), "w") as fh:
                fh.write(text.replace("{python}", sys.executable))
            os.chmod(os.path.join(b, name), 0o755)
        pids, self.beat = os.path.join(self.td, "pids"), os.path.join(self.td, "beat")
        env = dict(os.environ, PATH=b + os.pathsep + os.environ.get("PATH", ""),
                   FAKE_SKY_WRITES=os.pathsep.join(writes), FAKE_SKY_PIDS=pids,
                   FAKE_SKY_BEAT=self.beat, PYTHONDONTWRITEBYTECODE="1")
        for k in ("SLURM_JOB_PARTITION", "SLURM_JOB_ID"):
            env.pop(k, None)
        if limit:
            env.update(SLURM_JOB_ID="4242", FAKE_SCONTROL_LIMIT=limit)
        if kill_parent:
            env["FAKE_SKY_KILL_PARENT"] = "1"
        job = subprocess.Popen([sys.executable, os.path.join(SCRIPTS, "skyline_peaks.py")] + argv,
                               env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.pids.append(job.pid)
        deadline = time.time() + 60
        while not os.path.exists(pids):
            if job.poll() is not None:
                self.fail(f"the job ended before SkylineCmd started: {job.communicate()}")
            self.assertLess(time.time(), deadline, "SkylineCmd never started")
            time.sleep(0.05)
        with open(pids) as fh:
            child, grand = (int(x) for x in fh.read().split())
        self.pids += [child, grand]
        if send:
            os.kill(job.pid, send)
        _out, err = job.communicate(timeout=90)
        self.child, self.grand = child, grand
        return job.returncode, err

    def assert_nothing_runs_on(self):
        """SkylineCmd and the process it started are gone, and nothing writes any more."""
        import time
        for pid in (self.child, self.grand):
            deadline = time.time() + 10
            while True:
                try:
                    os.kill(pid, 0)
                except ProcessLookupError:
                    break
                self.assertLess(time.time(), deadline, f"process {pid} still runs")
                time.sleep(0.05)
        size = os.path.getsize(self.beat)
        time.sleep(0.5)
        self.assertEqual(os.path.getsize(self.beat), size, "something still writes")

    def join_setup(self):
        out = os.path.join(self.td, "output", "skyline")
        w = os.path.join(out, sp.WORK)
        for d in ("status", "logs"):
            os.makedirs(os.path.join(w, d))
        p = dict(plan_for(raw=["/d/r1.raw"], work=w), out=out, plan_id="p1")
        sp.write_json(os.path.join(w, "plan.json"), p)
        sp.write_json(sp.record_path(out), {"schema": sp.SCHEMA, "status": "planned",
                                            "plan_id": "p1"})
        with open(os.path.join(w, "jobs.txt"), "w") as fh:
            fh.write("build\t1\nimport\t2\njoin\t3\n")
        for name, body in (("build.json", {"ok": True, "plan_id": "p1"}),
                           ("import_1.json", {"rc": 0, "plan_id": "p1"})):
            sp.write_json(os.path.join(w, "status", name), body)
        self.assertEqual(sp.record_state(os.path.dirname(out))[0], "submitted")
        return out, w, p

    def check_join_stopped(self, out, rc, err, reason):
        self.assertEqual(rc, 1, err)
        rec = sp.load_record(out)
        self.assertEqual((rec["status"], rec["stage"]), ("skipped", "join"))
        self.assertEqual(rec["reason"], reason)
        self.assertNotIn("90 min", rec["reason"])          # never the script's default (H4)
        # no partial picture, table or .sky.zip beside it
        self.assertEqual([f for f in os.listdir(out) if sp.OWN_OUTPUTS.match(f)], [])
        # finalize --zip keeps only what scratch_files.prune keeps: the record alone
        files = [f for f in os.listdir(out) if os.path.isfile(os.path.join(out, f))]
        dirs = [d for d in os.listdir(out) if os.path.isdir(os.path.join(out, d))]
        self.assertEqual(scratch_files.prune(out, dirs, files)[0], [sp.RECORD])
        self.assertEqual(dirs, [])
        # and nothing says "submitted ... make the report again once they finish" any more
        output_dir = os.path.dirname(out)
        said = [sp.readme_rows(out, "output/skyline")[0][1],
                sp.manifest_line(output_dir, True)[2], sp.report_data(output_dir)["md"],
                sp.deliverable(output_dir)[2]]
        for text in said:
            self.assertNotIn(sp.STATE_WORDS["submitted"], text)
            self.assertNotIn("Not finished", text)
        for text in said[:2] + said[3:]:
            self.assertIn(rec["reason"], text)
        # the report: the reason as a client reads it (sp.client_reason, C2 biologist H2/L3)
        self.assertIn(sp.client_reason(rec["reason"]), said[2])
        self.assertTrue(sp.report_data(output_dir)["md"].startswith("**Not run.**"))
        self.assert_nothing_runs_on()
        with open(os.path.join(out, sp.WORK, "logs", "join.txt")) as fh:
            self.assertIn("ended by skyline_peaks.py: Stopped", fh.read())

    def test_a_join_warned_of_its_time_limit(self):
        import signal
        out, w, p = self.join_setup()
        writes = [os.path.join(out, "skyline_precursors.tsv"), os.path.join(out, "peaks_G1.png"),
                  os.path.join(out, "disagreements.tsv"), p["share_zip"]]
        rc, err = self.run_job(["join", "--work", w, "--plan-id", "p1"], writes,
                               send=signal.SIGUSR1, limit="00:03:00")
        # the limit Slurm applied (scontrol), not JOB_MINUTES' 90 (H4: --time=3 said 90 min)
        self.check_join_stopped(out, rc, err, "the join job ran out of its time limit (Slurm's "
                                              "--time, 3 min) before it finished")

    def test_a_join_sent_sigterm_mid_batch(self):
        import signal
        out, w, p = self.join_setup()
        writes = [os.path.join(out, "skyline_precursors.tsv"), p["share_zip"]]
        rc, err = self.run_job(["join", "--work", w, "--plan-id", "p1"], writes,
                               kill_parent=True)
        # no limit to read (not in a job): no number at all, never a default one
        self.check_join_stopped(out, rc, err, "Slurm stopped the join job before it finished: "
                                              "at its time limit (Slurm's --time), or by a cancel "
                                              "(scancel) or a preemption")

    def test_an_import_sent_sigterm_says_why_and_the_join_names_it(self):
        import signal
        out = os.path.join(self.td, "skyline")
        w = os.path.join(out, sp.WORK)
        for d in ("status", "logs"):
            os.makedirs(os.path.join(w, d))
        p = dict(plan_for(raw=["/d/r1.raw"], work=w), out=out, plan_id="p1")
        sp.write_json(os.path.join(w, "plan.json"), p)
        stem = os.path.splitext(p["doc"])[0]
        for f in (stem + ".sky", stem + sp.BUILT_TAG + ".sky"):
            with open(f, "w") as fh:
                fh.write("<built/>")
        part = os.path.join(os.path.dirname(sp.task_doc(p, 1)), "r1.raw.skyd")
        rc, err = self.run_job(["import", "--work", w, "--task", "1", "--plan-id", "p1"],
                               [part], kill_parent=True, limit="00:04:00")
        self.assertEqual(rc, 1, err)
        with open(os.path.join(w, "status", "import_1.json")) as fh:
            st = json.load(fh)
        self.assertEqual((st["rc"], st["reason"]),
                         (1, "Slurm stopped the import job before it finished: at its time limit "
                             "(Slurm's --time, 4 min), or by a cancel (scancel) or a preemption"))
        self.assertEqual(st["document"], stem + sp.BUILT_TAG + ".sky")
        self.assertFalse(os.path.exists(os.path.join(w, "r1.raw.skyd")))   # never joined
        self.assert_nothing_runs_on()
        # F3: the join's record names the run with that reason, not "ended without a result"
        sp.write_json(sp.record_path(out), {"schema": sp.SCHEMA, "status": "planned",
                                            "plan_id": "p1"})
        sp.write_json(os.path.join(w, "status", "build.json"), {"ok": True, "plan_id": "p1"})
        with mock.patch("sys.stderr", new=io.StringIO()):
            self.assertEqual(sp.main(["join", "--work", w, "--plan-id", "p1"]), 1)
        self.assertIn(f"r1.raw: {st['reason']}", sp.load_record(out)["reason"])

    def test_a_followup_warned_of_its_time_limit_leaves_no_partial_list(self):
        import signal
        out = os.path.join(self.td, "skyline")
        w = os.path.join(out, sp.WORK)
        dest = os.path.join(self.td, sp.FOLLOWUP_DIR)
        for d in (os.path.join(w, "prm"), os.path.join(w, "logs"), dest):
            os.makedirs(d)
        p = dict(plan_for(work=w), out=out, plan_id="p1")
        sp.write_json(os.path.join(w, "plan.json"), p)
        targets = sp.prm_plan({"family": "orbitrap"}, ["thermo_fusion", "thermo_srm"], True, 45)
        sp.write_json(os.path.join(w, "prm", "prm.json"), {"targets": targets, "dest": dest})
        with open(os.path.join(w, "prm", "prm_batch.txt"), "w") as fh:
            fh.write("--in=x\n")
        lists = [os.path.join(dest, t["file"]) for t in targets]
        rc, err = self.run_job(["prm-run", "--work", w, "--plan-id", "p1"], lists[:1],
                               send=signal.SIGUSR1, limit="1-02:03:00")
        self.assertEqual(rc, 1, err)
        self.assertEqual(sorted(os.listdir(dest)), ["README.md", "followup_assay.json"])
        why = ("the follow-up assay job ran out of its time limit (Slurm's --time, 1563 min) "
               "before it finished")
        with open(os.path.join(dest, "README.md")) as fh:
            readme = fh.read()
        self.assertIn(f"**Not made.** {why}.", readme)
        self.assertIn("no list here is a partial one", readme)
        with open(os.path.join(dest, "followup_assay.json")) as fh:
            rec = json.load(fh)
        self.assertEqual((rec["status"], rec["reason"]), ("not made", why))
        self.assert_nothing_runs_on()

    def test_a_build_stopped_records_why(self):
        import signal
        with tempfile.TemporaryDirectory() as td:
            w, p = TestImportTimeAndDocument().work(td)
            sp.write_json(sp.record_path(td), {"schema": sp.SCHEMA, "status": "planned",
                                               "plan_id": "p1"})
            why = sp.stopped_reason("build", signal.SIGTERM)

            def stopped(*a, **k):
                raise sp.Stopped(why)
            before = signal.getsignal(signal.SIGTERM)
            sel = {"mode": "de_hits", "targets": [{"protein_group": "P1", "gene": "G1"}]}
            with mock.patch.object(sp, "needs_for_pictures", lambda: None), \
                    mock.patch.object(sp, "verify_container", lambda plan: {"name": "Skyline"}), \
                    mock.patch.object(sp, "select_targets", lambda plan: dict(sel)), \
                    mock.patch.object(sp, "fasta_subset", lambda *a: {"P1"}), \
                    mock.patch.object(sp, "diann_rows", lambda *a: (4, 30.0)), \
                    mock.patch.object(sp, "run_skyline", stopped), \
                    mock.patch("sys.stderr", new=io.StringIO()):
                self.assertEqual(sp.main(["build", "--work", w, "--plan-id", "p1"]), 1)
            self.assertEqual(signal.getsignal(signal.SIGTERM), before)   # restored after the job
            rec = sp.load_record(td)
            self.assertEqual((rec["status"], rec["stage"], rec["reason"]), ("skipped", "build", why))
            with open(os.path.join(w, "status", "build.json")) as fh:
                self.assertEqual(json.load(fh), {"ok": False, "plan_id": "p1", "reason": why})

    def test_the_reasons(self):
        import signal
        with mock.patch.object(sp, "job_time_limit", return_value=None):
            self.assertEqual(sp.stopped_reason("import", signal.SIGTERM),
                             "Slurm stopped the import job before it finished: at its time limit "
                             "(Slurm's --time), or by a cancel (scancel) or a preemption")
        with mock.patch.dict(os.environ, {"SLURM_JOB_PARTITION": "low"}), \
                mock.patch.object(sp, "job_time_limit", return_value=(90 * 60, "scontrol")):
            low = sp.stopped_reason("join", signal.SIGUSR1)
        self.assertEqual(low, "the join job ran out of its time limit (Slurm's --time, 90 min) "
                              "before it finished, or Slurm preempted it; on the preemptible low "
                              "partition a preempted job is requeued and writes its own result")

    def scontrol(self, stdout, rc=0, job="77", which="/sw/slurm/bin/scontrol", exc=None):
        """job_time_limit() against a stand-in `scontrol show job -o`, and the calls made."""
        calls = []

        def run(argv, **k):
            calls.append((argv, k))
            if exc:
                raise exc
            return mock.Mock(returncode=rc, stdout=stdout, stderr="")
        env = {"SLURM_JOB_ID": job} if job else {}
        with mock.patch.dict(os.environ, env), \
                mock.patch.object(sp.shutil, "which", lambda name: which), \
                mock.patch.object(sp.subprocess, "run", run):
            if not job:
                os.environ.pop("SLURM_JOB_ID", None)
            return sp.job_time_limit(), calls

    def test_the_limit_slurm_applied_is_read(self):
        """H4 (Ruling 76): jobs resubmitted with --time=4 and --time=3 recorded "60 min" and
        "90 min", JOB_MINUTES' defaults. The limit is now Slurm's own, from scontrol."""
        import signal
        line = "JobId=77 JobName=sky_import RunTime=00:01:57 TimeLimit={} TimeMin=N/A\n"
        for shown, secs in (("00:04:00", 240), ("01:30:00", 5400), ("1-02:03:00", 93780),
                            ("00:04:30", 270)):
            with self.subTest(shown=shown):
                got, calls = self.scontrol(line.format(shown))
                self.assertEqual(got[0], secs)
                self.assertIn("scontrol", got[1])
                (argv, k), = calls
                self.assertEqual(argv, ["/sw/slurm/bin/scontrol", "show", "job", "-o", "77"])
                self.assertEqual(k.get("timeout"), sp.SCONTROL_TIMEOUT_S)
        with mock.patch.object(sp, "job_time_limit", return_value=(240, "scontrol")):
            self.assertEqual(sp.stopped_reason("import", signal.SIGUSR1),
                             "the import job ran out of its time limit (Slurm's --time, 4 min) "
                             "before it finished")
        with mock.patch.object(sp, "job_time_limit", return_value=(270, "scontrol")):
            self.assertIn("(Slurm's --time, 4 min 30 s)", sp.stopped_reason("import",
                                                                           signal.SIGUSR1))

    def test_no_number_when_the_limit_cannot_be_read(self):
        """Never a wrong number: when the limit cannot be read for certain, the reason says
        "its time limit" with none -- never JOB_MINUTES."""
        import signal
        import subprocess
        line = "JobId=77 RunTime=00:01:57 TimeLimit={} TimeMin=N/A\n"
        cases = {"not in a job": dict(stdout=line.format("00:04:00"), job=None),
                 "no scontrol": dict(stdout=line.format("00:04:00"), which=None),
                 "no answer in time": dict(stdout="", exc=subprocess.TimeoutExpired("scontrol", 10)),
                 "cannot start": dict(stdout="", exc=OSError("no such file")),
                 "an error": dict(stdout="slurm_load_jobs error: Invalid job id", rc=1),
                 "UNLIMITED": dict(stdout=line.format("UNLIMITED")),
                 "Partition_Limit": dict(stdout=line.format("Partition_Limit")),
                 "two limits": dict(stdout=line.format("00:04:00") + line.format("00:05:00")),
                 "no TimeLimit": dict(stdout="JobId=77 RunTime=00:01:57\n")}
        for name, kw in cases.items():
            with self.subTest(case=name):
                got, calls = self.scontrol(**kw)
                self.assertIsNone(got)
                if name in ("not in a job", "no scontrol"):
                    self.assertEqual(calls, [])
        with mock.patch.object(sp, "job_time_limit", return_value=None):
            for step in sp.JOB_MINUTES:
                for signum in (signal.SIGTERM, signal.SIGUSR1):
                    why = sp.stopped_reason(step, signum)
                    self.assertIn("its time limit (Slurm's --time)", why)
                    self.assertNotRegex(why, r"\d")

    def test_every_job_asks_for_the_early_warning_and_is_its_batch_shell(self):
        """--signal=B:USR1@120: "B:" reaches the batch shell only -- Slurm signals job steps
        without it, and these jobs have none -- and the batch shell IS skyline_peaks.py, because
        the job script exec's it."""
        p = plan_for()
        for queue in (("high", "a", "a-high-qos"), ("low", "publicgrp", "publicgrp-low-qos")):
            for step, sub in sp.JOB_SUBCOMMAND.items():
                with self.subTest(step=step, partition=queue[0]):
                    t = sp.job_script(p, step, queue, array="1-3%3" if step == "import" else None)
                    self.assertEqual(t.count("#SBATCH --signal=B:USR1@120"), 1)
                    self.assertEqual(sp.WARN_SECONDS, 120)
                    last = t.strip().splitlines()[-1]
                    self.assertTrue(last.startswith("exec "), last)
                    self.assertIn(f"skyline_peaks.py {sub} --work ", last)


class TestOneDocumentName(unittest.TestCase):
    """Rule 3: the Skyline document's name, <name>_top_hits.sky.zip, is defined once
    (skyline_peaks.SKYLINE_DOCUMENT_SUFFIX). core_submission had its own copy, and make_report and
    this script spelled it out again; each now reads the one definition."""

    @staticmethod
    def source(module):
        with open(os.path.splitext(module.__file__)[0] + ".py", encoding="utf-8") as fh:
            return fh.read()

    def test_core_submission_imports_it(self):
        import core_submission as cs
        # the same object, not an equal copy (a second literal is a second object)
        self.assertIs(cs.SKYLINE_DOCUMENT_SUFFIX, sp.SKYLINE_DOCUMENT_SUFFIX)
        self.assertNotIn('SKYLINE_DOCUMENT_SUFFIX = "', self.source(cs))
        self.assertEqual(cs.skyline_document({"skyline/peaks_G.png", "skyline/T_top_hits.sky.zip"}),
                         "skyline/T_top_hits.sky.zip")

    def test_no_other_copy_in_this_script_or_the_catalog(self):
        import make_report
        self.assertEqual(sp.SKYLINE_DOCUMENT_SUFFIX, "_top_hits.sky.zip")
        self.assertEqual(self.source(sp).count("_top_hits"), 1, "one definition in skyline_peaks")
        self.assertNotIn("_top_hits", self.source(make_report))
        self.assertEqual(make_report.describe("X" + sp.SKYLINE_DOCUMENT_SUFFIX)[0],
                         make_report.PEAKS)
        self.assertTrue(sp.OWN_OUTPUTS.match("X" + sp.SKYLINE_DOCUMENT_SUFFIX))
        self.assertFalse(sp.OWN_OUTPUTS.match("X_top_hitsXsky.zip"))     # an escaped pattern
        with tempfile.TemporaryDirectory() as td:
            search = TestNotRun().make_search(td)
            rc, out, said = plan_here(td, search)
            self.assertEqual(rc, 0, said)
            with open(os.path.join(out, sp.WORK, "plan.json")) as fh:
                plan = json.load(fh)
        self.assertEqual(os.path.basename(plan["share_zip"]), "T" + sp.SKYLINE_DOCUMENT_SUFFIX)
        self.assertEqual(os.path.basename(plan["doc"]) + ".zip", "T" + sp.SKYLINE_DOCUMENT_SUFFIX)

    def test_core_submission_spells_it_nowhere(self):
        """L9 review F9: core_submission still spelled `*_top_hits.sky.zip` in the email
        draft's warning (and `_top_hits` in a docstring and a comment); the warning now names
        the document by the one definition."""
        import core_submission as cs
        from test_core_submission import write_summary
        self.assertNotIn("_top_hits", self.source(cs))
        with tempfile.TemporaryDirectory() as td:
            summary, _s = write_summary(td)
            dj = os.path.join(td, "delivery.json")         # one from before 2.14.2: no `skyline`
            with open(dj, "w") as fh:
                json.dump({"contents": ["Analysis_Report.html", "skyline"], "mode": "analysis"},
                          fh)
            with mock.patch("sys.stdout", new=io.StringIO()) as so:
                cs.main(["email-draft", "--summary", summary, "--delivery", dj,
                         "--out", os.path.join(td, "EMAIL_DRAFT.md")])
            # the records lane's finalize-state warning also says "before skill 2.14.2" for
            # this delivery.json; this test is about the Skyline one only
            warned = [w for w in json.loads(so.getvalue())["warnings"]
                      if "before skill 2.14.2" in w and "Skyline" in w]
        self.assertEqual(len(warned), 1)
        self.assertIn(f"skyline/ folder for a *{sp.SKYLINE_DOCUMENT_SUFFIX} and add", warned[0])

    def test_the_docs_name_the_rule_and_what_status_prints(self):
        """L9 review F9: references/skyline.md named core_submission.skyline_deliverable as the
        delivery rule (since b8a04605 it is skyline_peaks.deliverable), and said `status` "says
        more" about a job that ran out (it prints the record and the job ids, not job states)."""
        with open(os.path.join(os.path.dirname(SCRIPTS), "references", "skyline.md"),
                  encoding="utf-8") as fh:
            doc = " ".join(fh.read().split())
        self.assertIn("`skyline_peaks.deliverable`, the one rule", doc)
        self.assertNotIn("say more", doc)
        self.assertIn("`status` lists its job ids", doc)
        with tempfile.TemporaryDirectory() as td:
            sky = os.path.join(td, "skyline")
            sp.write_json(sp.record_path(sky), {"schema": sp.SCHEMA, "status": "planned",
                                                "plan_id": "p", "written": "t"})
            os.makedirs(os.path.join(sky, sp.WORK))
            with open(os.path.join(sky, sp.WORK, "jobs.txt"), "w") as fh:
                fh.write("build\t11\njoin\t13\n")
            with mock.patch("sys.stdout", new=io.StringIO()) as so:
                sp.main(["status", "--out", sky])
            self.assertEqual(json.loads(so.getvalue())["jobs"], ["build", "11", "join", "13"])
            with self.assertRaises(SystemExit) as e:
                sp.main(["plan", "--search", td, "--out", sky])
            self.assertIn("`status` lists its job ids; `sacct -j <ids>` says whether they still "
                          "run", str(e.exception))


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
            # a finished check: only that is delivered (core_submission.skyline_deliverable)
            sp.write_json(sp.record_path(os.path.join(td, "skyline")), dict(OK_REC))
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
