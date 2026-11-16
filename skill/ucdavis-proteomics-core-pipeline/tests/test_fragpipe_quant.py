#!/usr/bin/env python3
"""The FragPipe DDA protein quantity (2.15 F4; Rulings 64, 90): fragpipe_quant.py.

  * the quantity columns are EXACTLY "<exp> MaxLFQ Intensity" for every experiment, else exactly
    "<exp> Intensity" -- never a mix, never a Unique / Total / Razor variant (2.14.2's adapter took
    any column ending in "Intensity");
  * a run of the map missing from combined_protein.tsv (or from directLFQ's table) stops, naming it;
  * 2 experiments -> directLFQ is the quantity, IonQuant Intensity beside it; IonQuant's own log
    line is quoted as why there is no MaxLFQ; no map -> a caveat;
  * the guard refuses pandas >= 3, a COPY flag True, and a 0 shift applied when one was found --
    with no size threshold (a real shift can be 0.021 log2) -- and compares tables BY COLUMN NAME;
  * the venv is found only by its builder's marker with these pins; without it the refusal names
    the builder; a record for another combined_ion.tsv is never reused;
  * one pin list: the builder reads it, and refuses a folder it did not finish.
Stdlib only (directLFQ itself runs in its own venv: the HIVE real-data check).
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

import fragpipe_quant as fq                     # noqa: E402

BASE = ["Protein", "Protein ID", "Entry Name", "Gene", "Combined Spectral Count",
        "Combined Total Spectral Count"]


def header(exps, maxlfq=True, unique_total=True):
    h = list(BASE)
    for e in exps:
        h += [f"{e} Spectral Count", f"{e} Unique Spectral Count", f"{e} Total Spectral Count",
              f"{e} Intensity"]
        if unique_total:
            h += [f"{e} Unique Intensity", f"{e} Total Intensity", f"{e} Razor Intensity"]
        if maxlfq:
            h += [f"{e} MaxLFQ Intensity"]
            if unique_total:
                h += [f"{e} Unique MaxLFQ Intensity", f"{e} MaxLFQ Total Intensity"]
    return h


def run_map(exps, runs=None):
    runs = runs or [e.replace("_", "-") for e in exps]
    return [{"experiment": e, "run": r, "input": f"/stage/{r}.d"} for e, r in zip(exps, runs)]


class ExactColumns(unittest.TestCase):
    def test_maxlfq_for_every_experiment(self):
        q = fq.quantity_columns(header(["A", "B", "C"]), ["A", "B", "C"])
        self.assertEqual(q["kind"], "MaxLFQ Intensity")
        self.assertEqual(q["columns"], {e: f"{e} MaxLFQ Intensity" for e in "ABC"})

    def test_unique_and_total_intensity_are_never_read(self):
        """2.14.2 took any column ending in 'Intensity': without MaxLFQ columns, 'A Unique
        Intensity' / 'A Total Intensity' became samples of their own."""
        h = header(["A", "B"], maxlfq=False)
        q = fq.quantity_columns(h, ["A", "B"])
        self.assertEqual(q["columns"], {"A": "A Intensity", "B": "B Intensity"})
        self.assertEqual(fq.experiments_from_header(h), ["A", "B"])
        for name in fq.experiments_from_header(header(["A", "B"])):
            self.assertNotIn("Unique", name)
            self.assertNotIn("Total", name)
        # an experiment whose only columns are the variants: stopped, never read
        h2 = [c for c in h if c != "B Intensity"]
        with self.assertRaises(fq.QuantError) as cm:
            fq.quantity_columns(h2, ["A", "B"])
        self.assertIn("B", str(cm.exception))

    def test_the_quantified_check_reads_only_quantity_columns(self):
        h = header(["A", "B"])
        self.assertEqual(sorted(fq.quantity_column_names(h)),
                         ["A Intensity", "A MaxLFQ Intensity", "B Intensity",
                          "B MaxLFQ Intensity"])

    def test_never_a_mix(self):
        h = [c for c in header(["A", "B", "C"], unique_total=False) if c != "C MaxLFQ Intensity"]
        q = fq.quantity_columns(h, ["A", "B", "C"])
        self.assertEqual(q["kind"], "Intensity")
        self.assertEqual(set(q["columns"].values()), {"A Intensity", "B Intensity", "C Intensity"})
        self.assertEqual(q["maxlfq_missing"], ["C"])

    def test_a_run_missing_from_the_table_stops_naming_it(self):
        h = header(["S1", "S2"])
        with self.assertRaises(fq.QuantError) as cm:
            fq.plan(h, run_map(["S1", "S2", "S3"]))
        self.assertIn("S3", str(cm.exception))
        self.assertIn("nothing was built", str(cm.exception))


class ThePlan(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def test_two_experiments_are_directlfq_with_ionquant_beside(self):
        with open(os.path.join(self.d, "log_2026-11-01_20-59-15.txt"), "w") as fh:
            fh.write("IonQuant ...\n[WARNING] - There are only 2 experiments. Will not calculate "
                     "MaxLFQ intensity.\n")
        p = fq.plan(header(["S_1", "S_2"], maxlfq=False), run_map(["S_1", "S_2"]), self.d)
        self.assertEqual(p["quantity"], fq.QUANTITY_DIRECTLFQ)
        self.assertEqual(p["n_experiments"], 2)
        self.assertEqual(p["companion_column"], fq.COMPANION_COLUMN)
        self.assertEqual(p["ionquant"]["columns"], {"S_1": "S_1 Intensity", "S_2": "S_2 Intensity"})
        self.assertIn("There are only 2 experiments. Will not calculate MaxLFQ intensity.",
                      p["why_not_maxlfq"])
        self.assertEqual([e["run"] for e in p["experiments"]], ["S-1", "S-2"])
        self.assertIsNone(p["caveat"])
        self.assertTrue(p["run_map"])

    def test_three_with_maxlfq_is_ionquant_maxlfq(self):
        p = fq.plan(header(["A", "B", "C"]), run_map(["A", "B", "C"]), self.d)
        self.assertEqual(p["quantity"], fq.QUANTITY_MAXLFQ)
        self.assertIsNone(p["why_not_maxlfq"])
        self.assertIsNone(p["companion_column"])

    def test_three_without_maxlfq_are_directlfq_and_say_why(self):
        """Ruling 102 (3): with no MaxLFQ columns, directLFQ at ANY run count -- IonQuant's
        Intensity is inflated by depth whatever the number of runs (Ruling 64) -- with the
        Intensity beside it, through the same declaration path."""
        with open(os.path.join(self.d, "fragpipe.workflow"), "w") as fh:
            fh.write("ionquant.maxlfq=0\n")
        for n in (3, 8):
            exps = [f"E{i}" for i in range(n)]
            p = fq.plan(header(exps, maxlfq=False), run_map(exps), self.d)
            self.assertEqual(p["quantity"], fq.QUANTITY_DIRECTLFQ, n)
            self.assertEqual(p["n_experiments"], n)
            self.assertEqual(p["companion_column"], fq.COMPANION_COLUMN)
            self.assertEqual(p["ionquant"]["columns"], {e: f"{e} Intensity" for e in exps})
            self.assertIn("ionquant.maxlfq=0", p["why_not_maxlfq"])
            self.assertIsNone(p["caveat"])

    def test_a_requested_directlfq_says_requested_not_two_experiments(self):
        """--directlfq (or an --overrides ionquant.maxlfq=0) with 3 or more experiments: the
        reason is the request, and the one phrase for the report says so."""
        with open(os.path.join(self.d, "fragpipe.workflow"), "w") as fh:
            fh.write("ionquant.maxlfq=0\n")
        exps = ["E1", "E2", "E3", "E4"]
        p = fq.plan(header(exps, maxlfq=False), run_map(exps), self.d)
        self.assertEqual(p["quantity"], fq.QUANTITY_DIRECTLFQ)
        self.assertTrue(p["why_not_maxlfq"].startswith("directLFQ was requested"))
        self.assertNotIn("only from 3", p["why_not_maxlfq"])
        self.assertEqual(fq.reason_phrase(p), "requested instead of IonQuant's MaxLFQ")

    def test_two_experiments_keep_the_true_two_experiment_reason(self):
        p = fq.plan(header(["S_1", "S_2"], maxlfq=False), run_map(["S_1", "S_2"]), self.d)
        self.assertNotIn("requested", p["why_not_maxlfq"])
        self.assertIn("2 experiments", fq.reason_phrase(p))
        with open(os.path.join(self.d, "fragpipe.workflow"), "w") as fh:
            fh.write("ionquant.maxlfq=0\n")
        p = fq.plan(header(["S_1", "S_2"], maxlfq=False), run_map(["S_1", "S_2"]), self.d)
        self.assertIn("requested", p["why_not_maxlfq"])
        self.assertIn("only from 3", p["why_not_maxlfq"])     # both are true; both said

    def test_some_maxlfq_columns_missing_is_directlfq_never_a_mix(self):
        h = [c for c in header(["A", "B", "C"]) if c != "C MaxLFQ Intensity"]
        p = fq.plan(h, run_map(["A", "B", "C"]), self.d)
        self.assertEqual(p["quantity"], fq.QUANTITY_DIRECTLFQ)
        self.assertIn("only some experiments", p["why_not_maxlfq"])

    def test_one_experiment_is_intensity_and_says_there_is_no_between_run_step(self):
        p = fq.plan(header(["A"], maxlfq=False), run_map(["A"]), self.d)
        self.assertEqual(p["quantity"], fq.QUANTITY_INTENSITY)
        self.assertIn("one experiment", p["why_not_maxlfq"])
        self.assertIn("no between-run", p["why_not_maxlfq"])
        self.assertIsNone(p["companion_column"])

    def test_no_map_is_a_caveat(self):
        p = fq.plan(header(["A", "B", "C"]), None, self.d)
        self.assertFalse(p["run_map"])
        self.assertEqual(p["caveat"], fq.NO_MAP_CAVEAT)
        self.assertEqual([e["experiment"] for e in p["experiments"]], ["A", "B", "C"])

    def test_a_map_naming_an_experiment_twice_is_refused(self):
        m = run_map(["A_1", "A_1", "B"], runs=["A-1", "A_1", "B"])
        with self.assertRaises(fq.QuantError) as cm:
            fq.plan(header(["A_1", "B"]), m, self.d)
        self.assertIn("A_1", str(cm.exception))

    def test_a_run_fragpipe_did_not_list_stops_naming_it(self):
        m = run_map(["A", "B", "C"])
        m[2]["experiment"] = None
        with self.assertRaises(fq.QuantError) as cm:
            fq.plan(header(["A", "B"]), m, self.d)
        self.assertIn("C", str(cm.exception))
        self.assertIn("did not search", str(cm.exception))

    def test_the_plan_holds_no_paths(self):
        p = fq.plan(header(["A", "B", "C"]), run_map(["A", "B", "C"]), self.d)
        self.assertNotIn("/stage/", json.dumps(p))


def table(rows):
    """{protein: {sample: intensity}} in read_protein_table's shape."""
    return {p: dict(v) for p, v in rows.items()}


PANDAS2 = {"directlfq": "0.3.3", "pandas": "2.2.3", "copy_flag": False}


class TheGuard(unittest.TestCase):
    """Ruling 90: pandas < 3 AND the COPY flag False AND the shift applied = the shift found."""

    def tables(self, shift_b, applied_b, n=40):
        """Unnormalised and normalised tables where sample B was moved by `applied_b` log2 (each
        protein also carries its own constant, as directLFQ's per-protein scaling does)."""
        raw, norm = {}, {}
        for i in range(n):
            a, b, c = 1e6 * (1 + i), 2e6 * (1 + i % 7), 1.5e6
            k = 1.0 + 0.01 * i                         # per-protein scale, same in every sample
            raw[f"P{i}"] = {"A": a, "B": b, "C": c}
            norm[f"P{i}"] = {"A": a * k, "B": b * k * 2 ** applied_b, "C": c * k}
        return {"A": 0.0, "B": shift_b, "C": 0.0}, norm, raw

    def test_applied_as_found_passes(self):
        found, norm, raw = self.tables(0.021, 0.021)
        g = fq.guard(PANDAS2, found, norm, raw)
        self.assertTrue(g["passed"], g)
        self.assertIsNone(g["refusal"])

    def test_a_small_shift_not_applied_is_refused(self):
        """C0: the real shift was only 0.021 log2, and pandas 3 dropped it: no size threshold."""
        found, norm, raw = self.tables(0.021, 0.0)
        g = fq.guard(PANDAS2, found, norm, raw)
        self.assertFalse(g["passed"])
        self.assertFalse(g["checks"]["shift_applied_equals_found"]["passed"])
        self.assertIn("was not applied", g["refusal"])
        self.assertIn("REFUSED", g["refusal"])

    def test_no_shift_found_and_none_applied_passes(self):
        found, norm, raw = self.tables(0.0, 0.0)
        self.assertTrue(fq.guard(PANDAS2, found, norm, raw)["passed"])

    def test_pandas_3_is_refused(self):
        found, norm, raw = self.tables(0.1, 0.1)
        g = fq.guard(dict(PANDAS2, pandas="3.1.1"), found, norm, raw)
        self.assertFalse(g["passed"])
        self.assertFalse(g["checks"]["pandas_below_3"]["passed"])
        self.assertIn("pandas 3.1.1", g["refusal"])

    def test_the_copy_flag_true_is_refused(self):
        found, norm, raw = self.tables(0.1, 0.1)
        g = fq.guard(dict(PANDAS2, copy_flag=True), found, norm, raw)
        self.assertFalse(g["passed"])
        self.assertFalse(g["checks"]["copy_flag_false"]["passed"])
        self.assertIn("COPY_NUMPY_ARRAYS_DERIVED_FROM_PANDAS = True", g["refusal"])
        g = fq.guard(dict(PANDAS2, copy_flag=None), found, norm, raw)    # not read = not passed
        self.assertFalse(g["passed"])

    def test_the_sentence_reads_each_check(self):
        """Review L3: guard_sentence() said "(< 3)" and "copy flag False" whatever the checks found;
        each claim now comes from its own check (a passed guard reads as before)."""
        found, norm, raw = self.tables(0.021, 0.021)

        def sentence(env, found=found, norm=norm):
            g = fq.guard(env, found, norm, raw)
            return fq.guard_sentence({"guard": g, "versions": {"pandas": env.get("pandas")}})

        ok = sentence(PANDAS2)
        self.assertTrue(ok.startswith(f"passed: pandas {PANDAS2['pandas']} (< 3); directLFQ's "
                                      f"copy flag False; its between-sample normalisation "
                                      f"applied as found"), ok)
        bad = sentence(dict(PANDAS2, pandas="3.1.1"))
        self.assertIn("pandas 3.1.1 (NOT < 3)", bad)
        self.assertNotIn("(< 3)", bad)
        bad = sentence(dict(PANDAS2, copy_flag=True))
        self.assertIn("copy flag NOT False", bad)
        self.assertNotIn("copy flag False", bad)
        bad = sentence(PANDAS2, found=self.tables(0.021, 0.0)[0], norm=self.tables(0.021, 0.0)[1])
        self.assertIn("normalisation NOT applied as found", bad)
        for b in (bad,):
            self.assertTrue(b.startswith("REFUSED"), b)
        # a record with no checks claims none of them
        none = fq.guard_sentence({"guard": {"passed": False}, "versions": {"pandas": "2.2.3"}})
        self.assertNotIn("(< 3)", none)
        self.assertNotIn("copy flag False", none)

    def test_tables_compared_by_column_name(self):
        """pandas 2 and 3 write directLFQ's columns in different orders (C0): a dict per protein
        keyed by sample name; read_protein_table keeps names whatever the order."""
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        a, b = os.path.join(d, "n.tsv"), os.path.join(d, "u.tsv")
        with open(a, "w") as fh:
            fh.write("protein\tA\tB\nP1\t100\t400\nP2\t10\t0\n")
        with open(b, "w") as fh:
            fh.write("protein\tB\tA\nP1\t200\t100\nP2\t0\t10\n")
        n, s = fq.read_protein_table(a)
        u, _ = fq.read_protein_table(b)
        self.assertEqual(s, ["A", "B"])
        self.assertEqual(n["P1"], {"A": 100.0, "B": 400.0})
        self.assertEqual(u["P1"], {"A": 100.0, "B": 200.0})
        self.assertEqual(n["P2"], {"A": 10.0})                    # 0 = not quantified
        g = fq.shift_check({"A": 0.0, "B": 1.0}, n, u)
        self.assertTrue(g["passed"], g)
        self.assertEqual(g["applied"]["B"], 1.0)

    def test_an_unpinned_directlfq_is_refused(self):
        found, norm, raw = self.tables(0.1, 0.1)
        self.assertFalse(fq.guard(dict(PANDAS2, directlfq="0.3.4"), found, norm, raw)["passed"])

    def test_different_samples_in_the_two_tables_are_refused(self):
        found, norm, raw = self.tables(0.1, 0.1)
        for v in raw.values():
            v.pop("C")
        self.assertFalse(fq.guard(PANDAS2, found, norm, raw)["passed"])

    def test_no_normaliser_shifts_are_refused(self):
        _, norm, raw = self.tables(0.1, 0.1)
        self.assertFalse(fq.guard(PANDAS2, {}, norm, raw)["passed"])


class TheVenv(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def make(self, pins=fq.DIRECTLFQ_PINS, marker=True):
        prefix = os.path.join(self.d, "venv")
        os.makedirs(os.path.join(prefix, "bin"))
        py = os.path.join(prefix, "bin", "python")
        os.symlink(sys.executable, py)
        if marker:
            with open(os.path.join(prefix, fq.VENV_MARKER), "w") as fh:
                json.dump({"pins": list(pins)}, fh)
        return prefix, py

    def test_found_by_its_marker(self):
        prefix, py = self.make()
        with mock.patch.dict(os.environ, {fq.VENV_ENV: prefix}):
            self.assertEqual(fq.venv_python()[0], py)

    def test_other_pins_or_no_marker_is_not_the_venv(self):
        for kw in ({"pins": ("directlfq==0.3.3", "pyyaml")}, {"marker": False}):
            shutil.rmtree(os.path.join(self.d, "venv"), ignore_errors=True)
            prefix, _ = self.make(**kw)
            with mock.patch.dict(os.environ, {fq.VENV_ENV: prefix}):
                py, why = fq.venv_python()
            self.assertIsNone(py)
            self.assertIn(prefix, why)

    def test_the_core_s_venv_is_named_in_one_constant(self):
        """Ruling 102: one versioned folder, named once (rule 3)."""
        self.assertEqual(fq.GROUP_VENV, "/quobyte/proteomics-grp/tools/directlfq/0.3.3-pandas2.2.3")
        hits = []
        for fn in sorted(os.listdir(SCRIPTS)):
            path = os.path.join(SCRIPTS, fn)
            if os.path.isfile(path):
                with open(path, errors="replace") as fh:
                    hits += [fn for line in fh if "tools/directlfq" in line]
        self.assertEqual(hits, ["fragpipe_quant.py"])
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop(fq.VENV_ENV, None)
            self.assertEqual(fq.venv_candidates()[0][0], fq.GROUP_VENV)

    def test_bigger_runs_get_a_bigger_node(self):
        small, big = fq.directlfq_resources(330_000, 2), fq.directlfq_resources(2_000_000_000, 40)
        self.assertEqual((small["cpus"], small["mem"]), (1, "8G"))
        self.assertEqual((big["cpus"], big["mem"], big["minutes"]), (4, "24G", 70))
        self.assertEqual(fq.directlfq_resources(1, 5000)["minutes"], 720)

    def test_a_set_env_var_never_falls_back(self):
        with mock.patch.dict(os.environ, {fq.VENV_ENV: os.path.join(self.d, "nowhere")}):
            self.assertEqual(len(fq.venv_candidates()), 1)


class RunDirectLfq(unittest.TestCase):
    """The adapter's side: a record for THIS combined_ion.tsv, or a run in the venv."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        with open(os.path.join(self.d, fq.ION_TABLE), "w") as fh:
            fh.write("Protein\tModified Sequence\tS_1 Intensity\tS_2 Intensity\n")

    def write_record(self, status="ok", passed=True, samples=("S_1", "S_2"), sha=None):
        os.makedirs(os.path.join(self.d, fq.DIRECTLFQ_DIR), exist_ok=True)
        out = os.path.join(fq.DIRECTLFQ_DIR, "combined_ion.tsv.protein_intensities.tsv")
        with open(os.path.join(self.d, out), "w") as fh:
            fh.write("protein\t" + "\t".join(samples) + "\nsp|P1|X\t" + "\t".join(
                "100" for _ in samples) + "\n")
        rec = {"status": status, "pins": list(fq.DIRECTLFQ_PINS), "samples": list(samples),
               "input_sha256": sha or fq.sha256_file(os.path.join(self.d, fq.ION_TABLE)),
               "outputs": {"normalised": out}, "versions": {"pandas": "2.2.3"},
               "guard": {"passed": passed, "refusal": None if passed else "REFUSED: test",
                         "checks": {}}}
        with open(os.path.join(self.d, fq.DIRECTLFQ_DIR, fq.DIRECTLFQ_RECORD), "w") as fh:
            json.dump(rec, fh)
        return rec

    def test_without_the_venv_the_refusal_names_the_builder(self):
        with mock.patch.dict(os.environ, {fq.VENV_ENV: os.path.join(self.d, "none")}):
            with self.assertRaises(fq.QuantError) as cm:
                fq.run_directlfq(self.d, ["S_1", "S_2"])
        self.assertIn(fq.BUILDER, str(cm.exception))
        self.assertIn("on_compute.py", str(cm.exception))

    def test_a_record_for_this_input_is_reused(self):
        self.write_record()
        runner = mock.Mock()
        rec = fq.run_directlfq(self.d, ["S_1", "S_2"], runner=runner)
        runner.assert_not_called()
        self.assertEqual(fq.directlfq_values(self.d, rec)[("sp|P1|X", "S_2")], 100.0)

    def test_a_record_for_another_input_is_never_reused(self):
        self.write_record(sha="0" * 64)
        with mock.patch.dict(os.environ, {fq.VENV_ENV: os.path.join(self.d, "none")}):
            with self.assertRaises(fq.QuantError):
                fq.run_directlfq(self.d, ["S_1", "S_2"], runner=mock.Mock())

    def venv(self):
        prefix = os.path.join(self.d, "venv")
        os.makedirs(os.path.join(prefix, "bin"))
        os.symlink(sys.executable, os.path.join(prefix, "bin", "python"))
        with open(os.path.join(prefix, fq.VENV_MARKER), "w") as fh:
            json.dump({"pins": list(fq.DIRECTLFQ_PINS)}, fh)
        return prefix

    def test_a_refused_run_stops_with_the_guard_s_words(self):
        prefix = self.venv()

        def runner(cmd, **kw):
            self.assertIn("on_compute.py", cmd[1])             # a login node -> srun
            self.assertIn("directlfq", cmd)
            self.write_record(status="refused", passed=False)
            return subprocess.CompletedProcess(cmd, 3, "", "refused")
        with mock.patch.dict(os.environ, {fq.VENV_ENV: prefix}):
            with self.assertRaises(fq.QuantError) as cm:
                fq.run_directlfq(self.d, ["S_1", "S_2"], runner=runner)
        self.assertIn("REFUSED: test", str(cm.exception))

    def test_a_refused_record_without_a_venv_says_the_refusal(self):
        self.write_record(status="refused", passed=False)
        with mock.patch.dict(os.environ, {fq.VENV_ENV: os.path.join(self.d, "none")}):
            with self.assertRaises(fq.QuantError) as cm:
                fq.run_directlfq(self.d, ["S_1", "S_2"])
        self.assertIn("REFUSED: test", str(cm.exception))
        self.assertIn(fq.BUILDER, str(cm.exception))

    def test_an_experiment_missing_from_directlfq_s_table_stops(self):
        self.write_record(samples=("S_1",))
        with self.assertRaises(fq.QuantError) as cm:
            fq.run_directlfq(self.d, ["S_1", "S_2"])
        self.assertIn("S_2", str(cm.exception))

    def test_no_ion_table_stops(self):
        os.remove(os.path.join(self.d, fq.ION_TABLE))
        with self.assertRaises(fq.QuantError) as cm:
            fq.run_directlfq(self.d, ["S_1", "S_2"])
        self.assertIn(fq.ION_TABLE, str(cm.exception))


class TheRecord(unittest.TestCase):
    def test_written_beside_the_search_and_into_its_provenance(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        with open(os.path.join(d, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "fragpipe"}, fh)
        p = fq.plan(header(["A", "B", "C"]), run_map(["A", "B", "C"]), d)
        fq.record(d, p)
        with open(os.path.join(d, fq.RECORD)) as fh:
            self.assertEqual(json.load(fh)["quantity"], fq.QUANTITY_MAXLFQ)
        with open(os.path.join(d, "search_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov[fq.PROVENANCE_KEY]["n_experiments"], 3)
        self.assertEqual(prov["engine"], "fragpipe")

    def test_a_refusal_is_recorded_without_the_core_s_paths(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        why = ("directLFQ's pinned venv is not installed (/quobyte/grp/tools/directlfq/v (the "
               "Core's group location): not there). Build it: python3 /home/u/skill/scripts/"
               "on_compute.py -- bash /home/u/skill/scripts/build_directlfq_venv.sh")
        fq.record(d, {"status": "refused", "refusal": why})
        with open(os.path.join(d, fq.RECORD)) as fh:
            text = json.load(fh)["refusal"]
        self.assertNotIn("/quobyte", text)
        self.assertNotIn("/home/u", text)
        self.assertIn("on_compute.py", text)
        self.assertIn("build_directlfq_venv.sh", text)


BUILDER = os.path.join(SCRIPTS, fq.BUILDER)            # the skill's own script, not a search job


class OnePinList(unittest.TestCase):

    def test_the_builder_reads_the_pins(self):
        with open(BUILDER) as fh:
            text = fh.read()
        self.assertNotIn("pandas==", text)
        self.assertNotIn("directlfq==", text)
        self.assertIn("fragpipe_quant.py\" pins", text)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fragpipe_quant.py"), "pins"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.stdout.split(), list(fq.DIRECTLFQ_PINS))
        self.assertIn("pandas==2.2.3", fq.DIRECTLFQ_PINS)

    @unittest.skipUnless(shutil.which("bash"), "bash")
    def test_the_builder_refuses_a_folder_it_did_not_finish(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        prefix = os.path.join(d, "half")
        os.makedirs(prefix)
        with open(os.path.join(prefix, "keep.txt"), "w") as fh:
            fh.write("x")
        env = dict(os.environ, PATH=os.pathsep.join([os.path.dirname(sys.executable),
                                                     "/usr/bin", "/bin"]))
        r = subprocess.run(["bash", BUILDER, "--prefix", prefix, "--python", sys.executable],
                           capture_output=True, text=True, timeout=120, env=env)
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("Nothing was removed", r.stderr)
        self.assertTrue(os.path.exists(os.path.join(prefix, "keep.txt")))


if __name__ == "__main__":
    unittest.main()
