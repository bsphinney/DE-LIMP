#!/usr/bin/env python3
"""
Orbitrap mass accuracy with no documented DIA-NN value is MEASURED, not extrapolated -- and a
level that HAS a documented value keeps it.

DIA-NN's README gives an Orbitrap resolution -> accuracy table for 240k, 120k, 60k and 30k
only. Real DIA methods sit outside it: both Orbitraps in the FRAN re-search pilot (2026-09-16)
acquire MS2 at 15,000 -- Exploris 480 at 60k/15k and 120k/15k, Fusion Lumos at 120k/15k.
estimate_params.py extrapolated a log-log fit of the table to 23.3 ppm, tagged it
EXTRAPOLATED, and pinned it for the whole cohort; it also stamped that EXTRAPOLATED tag on the
MS1 value, which was a documented tier (60k -> 10, 120k -> 7).

The README's item 6 of "Changing default settings" is its method for optimising these on the
data: run DIA-NN on several representative runs and review what it recommends. So:

  * estimate_params.py emits NO mass-accuracy flag for an Orbitrap with a level outside the
    table, and records the plan `measure_with_diann` in its rationale sidecar, with any level
    that does have a tier under `mass_accuracy_documented`. Neither flag is written on its own:
    DIA-NN 2.7.0 fixes BOTH levels when either is given ("automatic optimisation will not be
    performed as at least one of MS1/MS2 mass accuracies is user-provided"; HIVE srun 23528991,
    `--mass-acc-ms1 7` alone -> "Mass accuracy will be fixed to 2e-05 (MS2) and 7e-06 (MS1)"),
    so a lone `--mass-acc-ms1 7` would silently fix MS2 at 20 ppm.
  * The documented level is pinned as documented. Measuring it instead pinned MS1 4.2 at 120k,
    and every DIA-NN pass then logged "WARNING: the MS1 mass accuracy setting (4.2 ppm) deviates
    significantly from the value recommended (7 ppm) for the Orbitrap resolution of this run
    (120000)" (review of this branch, HIVE compare/probe_median).
  * The 5-step chain treats that -- and only that -- unpinned mass accuracy as recoverable.
  * The source text classify_instrument() hands resolve_defaults.py / make_presets.py names no
    route: those write it as the manifest's ppm_source for Radiant and FragPipe too.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import diann_parallel as dp  # noqa: E402
import estimate_params as ep  # noqa: E402
import run_search  # noqa: E402


def estimate_proc(d, instrument, ms1_res=None, ms2_res=None, overrides=None, name="params.cfg"):
    """Run estimate_params.py exactly as SKILL.md step 6b does. Returns (process, cfg path)."""
    cfg = os.path.join(d, name)
    argv = [sys.executable, os.path.join(SCRIPTS, "estimate_params.py"), "--engine", "diann",
            "--acquisition", "DIA", "--instrument", instrument,
            "--precursor-mz-range", "350", "1201", "--out", cfg]
    if ms1_res:
        argv += ["--ms1-resolution", str(ms1_res)]
    if ms2_res:
        argv += ["--ms2-resolution", str(ms2_res)]
    if overrides:
        argv += ["--overrides", json.dumps(overrides)]
    return subprocess.run(argv, capture_output=True, text=True, timeout=60), cfg


def estimate(d, instrument, ms1_res=None, ms2_res=None, overrides=None, name="params.cfg"):
    """estimate_proc(), which must succeed. Returns (cfg path, sidecar)."""
    p, cfg = estimate_proc(d, instrument, ms1_res, ms2_res, overrides, name)
    if p.returncode != 0:
        raise AssertionError(p.stderr)
    with open(cfg + ".rationale.json") as fh:
        return cfg, json.load(fh)


def flags_in(cfg):
    return [ln.split()[0] for ln in open(cfg) if ln.strip()]


class Args:
    """The subset of the argparse namespace parallel_decision reads."""
    parallel_threshold = 5
    no_parallel = False


FILES = ["f%d.raw" % i for i in range(310)]


def _inputs(d, n=6):
    """Real input paths for the generator: .mzML files, so no .NET 8 install is attempted."""
    raws = []
    for i in range(n):
        raws.append(os.path.join(d, "run%d.mzML" % i))
        open(raws[-1], "w").close()
    fasta = os.path.join(d, "db.fasta")
    with open(fasta, "w") as fh:
        fh.write(">sp|P1|X\nPEPTIDER\n")
    return raws, fasta


def _generate(d, cfg, raws, fasta, *extra):
    """diann_parallel.py exactly as run_search.run_diann_parallel runs it: a subprocess."""
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                           "--diann", "/bin/true", "--raw", *raws, "--fasta", fasta,
                           "--out", os.path.join(d, "out"), "--cfg", cfg, *extra],
                          capture_output=True, text=True, timeout=120)


class EstimateParamsTests(unittest.TestCase):

    def test_ms2_at_15k_is_measured_and_the_documented_ms1_tier_is_kept(self):
        """The pilot's Exploris 480 (120k/15k) and Lumos (120k/15k) methods."""
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "Orbitrap Exploris 480", 120000, 15000)
            flags = flags_in(cfg)
            self.assertNotIn("--mass-acc", flags, "a 15k MS2 tolerance was pinned from no evidence")
            self.assertNotIn("--mass-acc-ms1", flags,
                             "a lone --mass-acc-ms1 makes DIA-NN fix MS2 at 20 ppm as well")
            self.assertEqual(side["mass_accuracy_plan"], ep.MEASURE_WITH_DIANN)
            self.assertEqual(side["mass_accuracy_documented"], {"--mass-acc-ms1": 7},
                             "120k MS1 has a README tier; it must not be measured away")
            text = json.dumps(side)
            self.assertNotIn("EXTRAPOLATED", text)
            self.assertNotIn("23.3", text)
            ms1 = side["rationale"]["--mass-acc-ms1"]
            self.assertEqual(ms1["value"], 7)
            self.assertIn("120,000", ms1["source"])
            self.assertIn(ep.SRC_TABLE, ms1["source"])
            ms2 = side["rationale"]["--mass-acc"]["source"]
            self.assertIn("15,000", ms2, "the source must say which level has no tier")
            self.assertIn("measured", ms2)

    def test_ms1_at_60k_ms2_at_15k_keeps_ms1_10(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "Orbitrap Exploris 480", 60000, 15000)
            self.assertNotIn("--mass-acc", flags_in(cfg))
            self.assertNotIn("--mass-acc-ms1", flags_in(cfg))
            self.assertEqual(side["mass_accuracy_plan"], ep.MEASURE_WITH_DIANN)
            self.assertEqual(side["mass_accuracy_documented"], {"--mass-acc-ms1": 10})

    def test_both_levels_outside_the_table_are_both_measured(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "Orbitrap Exploris 480", 15000, 15000)
            self.assertEqual(side["mass_accuracy_plan"], ep.MEASURE_WITH_DIANN)
            self.assertEqual(side["mass_accuracy_documented"], {})

    def test_classify_instrument_keeps_the_documented_level_and_names_no_route(self):
        """resolve_defaults.py and make_presets.py write this source into the manifest for EVERY
        engine; SKILL.md tells the agent to quote it. It said 'measured in step 1b of the 5-step
        chain' for Radiant, FragPipe and single-shot DIA-NN, none of which have a step 1b."""
        cls, ms1, ms2, _, src = ep.classify_instrument("Orbitrap Exploris 480", 120000, 15000)
        self.assertEqual((cls, ms1, ms2), ("orbitrap_untabled", 7, None))
        for s in (src, ep.classify_instrument("Orbitrap Fusion Lumos")[4]):
            self.assertNotIn("step 1b", s)
            self.assertNotIn("5-step", s)
        self.assertIn("15,000", src)
        self.assertIn("120,000", src)

    def test_ppm_for_resolution_does_not_extrapolate_outside_the_table(self):
        for res in (15000, 7500, 480000):
            ppm, src = ep.ppm_for_resolution(res)
            self.assertIsNone(ppm, f"{res}: extrapolated to {ppm}")
            self.assertIn("30k-240k", src)
        self.assertEqual(ep.ppm_for_resolution(120000)[0], 7)
        self.assertEqual(ep.ppm_for_resolution(30000)[0], 15)

    def test_documented_tiers_are_pinned_and_never_tagged_extrapolated(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "Orbitrap Exploris 480", 120000, 30000)
            lines = open(cfg).read().splitlines()
            self.assertIn("--mass-acc 15", lines)
            self.assertIn("--mass-acc-ms1 7", lines)
            self.assertEqual(side["mass_accuracy_plan"], "pinned")
            self.assertNotIn("EXTRAPOL", json.dumps(side))
            self.assertNotIn("interpolated", json.dumps(side))

    def test_each_level_carries_its_own_source(self):
        """The old code tagged BOTH flags with one level's source: MS2's when it was
        extrapolated, otherwise MS1's -- so a documented MS1 read EXTRAPOLATED, and an
        interpolated MS2 read as documented."""
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "Orbitrap Exploris 480", 120000, 45000)
            ms2_src = side["rationale"]["--mass-acc"]["source"]
            ms1_src = side["rationale"]["--mass-acc-ms1"]["source"]
            self.assertIn("interpolated", ms2_src)
            self.assertNotIn("interpolated", ms1_src, ms1_src)
            self.assertIn("120,000", ms1_src)

    def test_orbitrap_with_unknown_resolution_is_not_measured_because_ms1_would_be(self):
        """An Orbitrap whose RESOLUTION is unknown -- the default Thermo path, since a .raw
        carries no resolution here and a Thermo mzML usually has no MS:1000800 -- must NOT plan a
        measurement. With no resolution neither level has a tier, so BOTH would be measured, and
        this branch's own evidence says a measured MS1 is wrong: pinning the measured 4.2 ppm at
        120k made every DIA-NN pass log "the MS1 mass accuracy setting (4.2 ppm) deviates
        significantly from the value recommended (7 ppm) ... (120000)". It falls to auto instead,
        which the 5-step chain declines out loud."""
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "Orbitrap Fusion Lumos")
            self.assertNotIn("--mass-acc", flags_in(cfg))
            self.assertNotIn("--mass-acc-ms1", flags_in(cfg))
            self.assertEqual(side["mass_accuracy_plan"], ep.PLAN_AUTO)
            self.assertNotIn("orbitrap_generic", ep.MEASURE_CLASSES)
            self.assertIsNone(dp.mass_acc_measure_plan(cfg),
                              "an unknown-resolution Orbitrap was accepted for measurement")
            safe = dp.parallel_safe(cfg)
            self.assertFalse(safe["ok"])
            self.assertIn("--ms1-resolution", safe["remedy"] + safe["reason"])
            # ...and a resolution turns it straight back into a measurement with a documented MS1
            cfg2, side2 = estimate(d, "Orbitrap Fusion Lumos", 120000, 15000, name="res.cfg")
            self.assertEqual(side2["mass_accuracy_plan"], ep.MEASURE_WITH_DIANN)
            self.assertEqual(side2["mass_accuracy_documented"], {"--mass-acc-ms1": 7})

    def test_other_instruments_are_unchanged(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "timsTOF Pro", name="tims.cfg")
            self.assertIn("--mass-acc 15", open(cfg).read().splitlines())
            self.assertEqual(side["mass_accuracy_plan"], "pinned")
            cfg, side = estimate(d, "Astral", name="astral.cfg")
            self.assertIn("--mass-acc 10", open(cfg).read().splitlines())
            # not an Orbitrap: nothing documented and nothing we know how to measure for it
            cfg, side = estimate(d, "Mystery-9000", name="unknown.cfg")
            self.assertNotIn("--mass-acc", flags_in(cfg))
            self.assertEqual(side["mass_accuracy_plan"], "auto")

    def test_a_validated_sop_override_pins_it(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "Orbitrap Exploris 480", 120000, 15000,
                                 overrides={"--mass-acc": 20, "--mass-acc-ms1": 7})
            lines = open(cfg).read().splitlines()
            self.assertIn("--mass-acc 20", lines)
            self.assertEqual(side["mass_accuracy_plan"], "pinned")

    def test_a_one_flag_override_writes_the_other_levels_table_value_too(self):
        """references/parameters.md's own override example, `{"--mass-acc": 8}`, at the pilot's
        120k/15k wrote ONLY `--mass-acc 8`, while the sidecar still planned measure_with_diann
        with MS1 7 documented. DIA-NN 2.7.0 fixes the other level at 20 ppm when one flag is given
        (HIVE srun 23528991), the chain declined the cfg as mass_acc_unset, and a single-shot
        search ran MS1 at 20 ppm with no probe. origin/main wrote `--mass-acc-ms1 7 --mass-acc 8`
        for the same input. The overridden level is given and the other has a table value, so
        both are known: both are written, and nothing is left to measure."""
        cases = (((120000, 15000), {"--mass-acc": 8}, {"--mass-acc": "8", "--mass-acc-ms1": "7"}),
                 ((15000, 120000), {"--mass-acc-ms1": 9},
                  {"--mass-acc": "7", "--mass-acc-ms1": "9"}),
                 ((90000, 15000), {"--mass-acc": 8},
                  {"--mass-acc": "8", "--mass-acc-ms1": "7.5"}),
                 ((120000, 60000), {"--mass-acc": 8}, {"--mass-acc": "8", "--mass-acc-ms1": "7"}))
        for (r1, r2), ov, want in cases:
            with self.subTest(res=(r1, r2), overrides=ov), tempfile.TemporaryDirectory() as d:
                cfg, side = estimate(d, "Orbitrap Exploris 480", r1, r2, overrides=ov)
                lines = [ln.split(None, 1) for ln in open(cfg).read().splitlines()
                         if ln.split()[:1] in (["--mass-acc"], ["--mass-acc-ms1"])]
                self.assertEqual(len(lines), 2, lines)
                self.assertEqual(dict(lines), want)
                self.assertEqual(side["mass_accuracy_plan"], "pinned")
                self.assertEqual(side["mass_accuracy_documented"], {})
                (other,) = [f for f in want if f not in ov]
                src = side["rationale"][other]["source"]
                self.assertIn(ep.SRC_TABLE, src)
                self.assertNotIn(ep.LONE_FLAG_NOTE, src)
                for f in ov:
                    self.assertEqual(side["rationale"][f]["source"],
                                     "user-override (validated SOP)")
                self.assertIsNone(dp.mass_acc_measure_plan(cfg))
                safe = dp.parallel_safe(cfg)
                self.assertTrue(safe["ok"], safe["reason"])
                self.assertEqual(safe["measure"], ["window"])

    def test_a_one_flag_override_with_no_value_for_the_other_level_is_refused(self):
        """With no value for the other level -- outside the table, resolution unknown, or an
        unidentified instrument -- a one-flag override can only be written as a lone flag, and
        DIA-NN 2.7.0 then fixes the other level at 20 ppm. origin/main already did that silently
        for an Orbitrap of unknown resolution and for an unidentified instrument. Nothing is
        written; the message names the missing flag and what to give instead."""
        cases = (("MS1 overridden, MS2 15k", ("Orbitrap Exploris 480", 120000, 15000),
                  {"--mass-acc-ms1": 5}, "--mass-acc ", "measure"),
                 ("resolution unknown", ("Orbitrap Fusion Lumos", None, None),
                  {"--mass-acc": 8}, "--mass-acc-ms1", "--ms1-resolution"),
                 ("unidentified instrument", ("Mystery-9000", None, None),
                  {"--mass-acc": 8}, "--mass-acc-ms1", None))
        for label, (instr, r1, r2), ov, missing, remedy in cases:
            with self.subTest(label), tempfile.TemporaryDirectory() as d:
                p, cfg = estimate_proc(d, instr, r1, r2, overrides=ov)
                self.assertNotEqual(p.returncode, 0, f"{label}: {open(cfg).read()}"
                                    if os.path.exists(cfg) else label)
                self.assertNotIn("Traceback", p.stderr)
                self.assertIn(missing, p.stderr)
                self.assertIn("20 ppm", p.stderr)
                if remedy:
                    self.assertIn(remedy, p.stderr)
                self.assertFalse(os.path.exists(cfg), "a lone-flag cfg was written")
                self.assertFalse(os.path.exists(cfg + ".rationale.json"))
        with tempfile.TemporaryDirectory() as d:     # an EARLIER cfg at --out must not survive
            # "no cfg was written" is true of this run -- build_diann raises before the write --
            # but the previous run of the same command wrote one, and every later step reads the
            # path, not the mtime. A stale cfg beside a refusal on stdout, still carrying a
            # mass_accuracy_plan, is the worst of both: the chain would measure and pin for a
            # parameter set nobody approved.
            cfg, side = estimate(d, "Orbitrap Exploris 480", 120000, 30000)
            self.assertEqual(side["mass_accuracy_plan"], "pinned")
            self.assertTrue(os.path.exists(cfg) and os.path.exists(cfg + ".rationale.json"))
            p, cfg2 = estimate_proc(d, "Orbitrap Exploris 480", 120000, 15000,
                                    overrides={"--mass-acc-ms1": 5})
            self.assertEqual(cfg2, cfg)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("no cfg was written", p.stderr)
            self.assertFalse(os.path.exists(cfg), "an earlier cfg survived a refusal")
            self.assertFalse(os.path.exists(cfg + ".rationale.json"),
                             "an earlier rationale sidecar survived a refusal")

        with tempfile.TemporaryDirectory() as d:              # both given: nothing is missing
            cfg, side = estimate(d, "Orbitrap Fusion Lumos",
                                 overrides={"--mass-acc": 8, "--mass-acc-ms1": 5})
            lines = open(cfg).read().splitlines()
            self.assertIn("--mass-acc 8", lines)
            self.assertIn("--mass-acc-ms1 5", lines)
            self.assertEqual(side["mass_accuracy_plan"], "pinned")
            self.assertEqual(dp.parallel_safe(cfg)["measure"], ["window"])


class OtherRoutesTests(unittest.TestCase):
    """resolve_defaults.py -> workflow.manifest.json, and make_presets.py for Radiant, at the
    pilot's 120k/15k. Before this branch: ms1_ppm 7, ms2_ppm 23.3 (extrapolated), Radiant MS1
    extraction width 7. The first cut of this branch recorded ms1_ppm/ms2_ppm None with a
    ppm_source promising step 1b for every engine, and Radiant fell back to vendor 20/20."""

    def _resolve(self, d, engine, *more, name=None):
        dest = os.path.join(d, name or engine)
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                            "--engine", engine, "--acquisition", "DIA",
                            "--instrument", "Orbitrap Exploris 480",
                            "--ms1-res", "120000", "--ms2-res", "15000", "--dest", dest, *more],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(dest, "workflow.manifest.json")) as fh:
            return dest, json.load(fh)["search"]

    def test_a_lone_sop_ppm_flag_does_not_leave_its_level_described_as_measured(self):
        """ppm_source was gated on --ms1-ppm alone, so a site SOP that set only --ms2-ppm put the
        SOP's number in ms2_ppm while ppm_source still said MS2 "is measured with DIA-NN before
        the search". SKILL.md tells the agent to read that line out to the user for confirmation,
        so the user would have confirmed a measurement that was never going to happen."""
        with tempfile.TemporaryDirectory() as d:
            _, s = self._resolve(d, "diann", "--ms2-ppm", "20", name="ms2only")
            self.assertEqual((s["ms1_ppm"], s["ms2_ppm"]), (7, 20.0))
            self.assertIn("MS2: site SOP override", s["ppm_source"])
            self.assertNotIn("measured with DIA-NN", s["ppm_source"],
                             "the SOP supplied the only level there was to measure")
            # the MS1 half is still described, and it is the one from the table
            self.assertIn("120,000", s["ppm_source"])
            # a lone --ms1-ppm is the mirror image: MS2 is still the measured level
            _, s = self._resolve(d, "diann", "--ms1-ppm", "9", name="ms1only")
            self.assertEqual((s["ms1_ppm"], s["ms2_ppm"]), (9.0, None))
            self.assertIn("MS1: site SOP override", s["ppm_source"])
            self.assertIn("measured with DIA-NN", s["ppm_source"])
            # both given: nothing is derived at all
            _, s = self._resolve(d, "diann", "--ms1-ppm", "9", "--ms2-ppm", "20", name="both")
            self.assertEqual(s["ppm_source"], "site SOP override")
            # and with neither, nothing changes
            _, s = self._resolve(d, "diann", name="neither")
            self.assertIn("measured with DIA-NN", s["ppm_source"])
            self.assertNotIn("SOP", s["ppm_source"])

    def test_the_manifest_keeps_ms1_and_names_what_each_engine_does_with_ms2(self):
        with tempfile.TemporaryDirectory() as d:
            for engine in ("diann", "radiant", "fragpipe"):
                _, s = self._resolve(d, engine)
                self.assertEqual((s["ms1_ppm"], s["ms2_ppm"]), (7, None), engine)
                self.assertIn("15,000", s["ppm_source"], engine)
                if engine != "diann":
                    self.assertNotIn("step 1b", s["ppm_source"], engine)
                    self.assertNotIn("measured", s["ppm_source"], engine)
            _, s = self._resolve(d, "diann")
            self.assertIn("measured with DIA-NN", s["ppm_source"])

    def test_radiant_keeps_the_documented_ms1_width_and_says_ms2_is_vendor(self):
        with tempfile.TemporaryDirectory() as d:
            dest, s = self._resolve(d, "radiant")
            text = open(s["params_file"]).read()
            self.assertRegex(text, r"ms1ExtractionWidthPPM = 7(\.0)?\n")
            self.assertRegex(text, r"ms2ExtractionWidthPPM = 20(\.0)?\n")
            prov = s["preset_provenance"]
            self.assertIn("MS2", prov["tolerances"])
            self.assertIn("vendor", prov["tolerances"])
            self.assertNotIn("step 1b", json.dumps(prov))


class RoutingTests(unittest.TestCase):
    """Unpinned mass accuracy is recoverable ONLY when step 1b will measure it: the cfg omits
    both flags, its estimate_params.py sidecar says `measure_with_diann`, and the chain has a
    step 1b (probing on, no seed library)."""

    def setUp(self):
        self._real = run_search.slurm_available
        run_search.slurm_available = lambda: True

    def tearDown(self):
        run_search.slurm_available = self._real

    def test_measure_cfg_routes_to_the_chain(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "Orbitrap Exploris 480", 120000, 15000)
            use, why = run_search.parallel_decision("diann", FILES, cfg, Args())
            self.assertTrue(use, why)
            self.assertIn("mass accuracy", why)
            self.assertIn("step 1b", why)
            safe = dp.parallel_safe(cfg)
            self.assertTrue(safe["probe"])
            self.assertEqual(sorted(safe["measure"]), ["mass-acc", "window"])
            self.assertEqual(safe["mass_acc_documented"], {"--mass-acc-ms1": 7})

    def test_the_same_flags_without_the_plan_still_decline(self):
        """A hand-written cfg that simply forgot mass accuracy is a mistake to report -- on a
        timsTOF the documented 15/15 is right, and measuring would quietly replace it. The fix
        named says where the plan lives, so a cfg copied without its sidecar is explained."""
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "Orbitrap Exploris 480", 120000, 15000)
            os.remove(cfg + ".rationale.json")
            use, why = run_search.parallel_decision("diann", FILES, cfg, Args())
            self.assertFalse(use)
            self.assertIn("mass accuracy is not pinned", why)
            self.assertEqual(dp.parallel_safe(cfg)["code"], "mass_acc_unset")
            self.assertIn(".rationale.json", why)

    def test_unknown_instrument_still_declines(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "Mystery-9000")
            use, why = run_search.parallel_decision("diann", FILES, cfg, Args())
            self.assertFalse(use)

    def test_the_plan_never_rescues_a_zero_partial_or_junk_value(self):
        """Upstream's order holds: an invalid value is reported as invalid (never as a plan to
        measure, and never as the override-able `mass_acc_unset`), and a partial pin -- which
        DIA-NN 2.7.0 turns into BOTH levels fixed -- is not a plan either."""
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "Orbitrap Exploris 480", 120000, 15000)
            base = open(cfg).read()
            for label, extra, code in (
                    ("zero", "--mass-acc 0\n--mass-acc-ms1 0\n", "mass_acc_invalid"),
                    ("ms1 only", "--mass-acc-ms1 10\n", "mass_acc_unset"),
                    ("ms2 only", "--mass-acc 20\n", "mass_acc_unset"),
                    ("junk", "--mass-acc wide\n--mass-acc-ms1 7\n", "mass_acc_invalid"),
                    ("negative", "--mass-acc -3\n--mass-acc-ms1 7\n", "mass_acc_invalid"),
                    ("window junk", "--window wide\n", "window_invalid"),
                    ("window fractional", "--window 0.5\n", "window_invalid")):
                with open(cfg, "w") as fh:
                    fh.write(base + extra)
                safe = dp.parallel_safe(cfg)
                self.assertFalse(safe["ok"], label)
                self.assertEqual(safe["code"], code, label)
                if code != "window_invalid":        # there the plan stands; the window is junk
                    self.assertIsNone(dp.mass_acc_measure_plan(cfg), label)
                use, _ = run_search.parallel_decision("diann", FILES, cfg, Args())
                self.assertFalse(use, label)

    def test_a_sidecar_whose_documented_value_is_not_a_positive_number_declines(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "Orbitrap Exploris 480", 120000, 15000)
            for junk in (0, -7, "seven", None, [7]):
                side["mass_accuracy_documented"] = {"--mass-acc-ms1": junk}
                with open(cfg + ".rationale.json", "w") as fh:
                    json.dump(side, fh)
                self.assertFalse(dp.parallel_safe(cfg)["ok"], junk)
                self.assertIsNone(dp.mass_acc_measure_plan(cfg), junk)

    def test_no_step1b_means_nothing_measures_it(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "Orbitrap Exploris 480", 120000, 15000)
            s = dp.parallel_safe(cfg, seed_lib=os.path.join(d, "seed.parquet"))
            self.assertFalse(s["ok"])
            self.assertEqual(s["code"], "mass_acc_seeded")
            self.assertIn("mass accuracy", s["reason"])
            self.assertIn("--seed-lib", s["remedy"])
            s = dp.parallel_safe(cfg, probe_window=False)
            self.assertFalse(s["ok"])
            self.assertEqual(s["code"], "mass_acc_no_probe")
            self.assertIn("mass accuracy", s["reason"])
            self.assertIn("--no-probe-window", s["remedy"])

    def test_the_override_covers_omitted_mass_accuracy_with_no_step1b_and_nothing_else(self):
        """--allow-auto-mass-acc overrides only OMITTED mass accuracy (upstream #70). A planned
        cfg the chain cannot measure (--no-probe-window) IS omitted mass accuracy, so the override
        applies, with its warning; it never applies to an invalid value in the same cfg."""
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "Orbitrap Exploris 480", 120000, 15000)
            raws, fasta = _inputs(d)
            p = _generate(d, cfg, raws, fasta, "--no-probe-window")
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("--no-probe-window", p.stderr)
            self.assertIn("--allow-auto-mass-acc", p.stderr)
            p = _generate(d, cfg, raws, fasta, "--no-probe-window", "--allow-auto-mass-acc")
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("WARNING", p.stderr)
            info = json.loads(p.stdout)
            self.assertEqual(info["step1b_measures"], [])
            self.assertFalse(info["mass_acc"]["fixed"])
            with open(cfg, "a") as fh:
                fh.write("--mass-acc 0\n--mass-acc-ms1 0\n")
            p = _generate(d, cfg, raws, fasta, "--allow-auto-mass-acc")
            self.assertNotEqual(p.returncode, 0, "an invalid mass accuracy was overridden")
            self.assertIn("0 ppm", p.stderr)

    def test_a_pinned_window_still_gets_step1b_for_mass_accuracy(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "Orbitrap Exploris 480", 120000, 15000)
            with open(cfg, "a") as fh:
                fh.write("--window 7\n")
            s = dp.parallel_safe(cfg)
            self.assertTrue(s["ok"], s["reason"])
            self.assertTrue(s["probe"])
            self.assertEqual(s["measure"], ["mass-acc"])

    def test_an_interpolated_level_is_not_called_documented(self):
        """MS1 at 90k is between tiers, so its 7.5 ppm is interpolated from DIA-NN's table -- not a
        README value. The gate's reason and the chain's provenance said '--mass-acc-ms1 7.5 as
        documented'; only the level's own rationale entry said 'interpolated'."""
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "Orbitrap Exploris 480", 90000, 15000)
            self.assertEqual(side["mass_accuracy_documented"], {"--mass-acc-ms1": 7.5})
            self.assertIn("interpolated for 90,000",
                          side["rationale"]["mass_accuracy_documented"]["source"])
            cfg120, side120 = estimate(d, "Orbitrap Exploris 480", 120000, 15000, name="t.cfg")
            self.assertNotIn("interpolated", side120["rationale"]["mass_accuracy_documented"]
                             ["source"], "an exact 120k tier was described as interpolated")
            safe = dp.parallel_safe(cfg)
            self.assertIn("7.5", safe["reason"])
            self.assertNotIn("as documented", safe["reason"])
            self.assertIn("interpolated", safe["reason"])
            raws, fasta = _inputs(d)
            p = _generate(d, cfg, raws, fasta)
            self.assertEqual(p.returncode, 0, p.stderr)
            src = json.loads(p.stdout)["mass_acc"]["source"]
            self.assertIn("7.5", src)
            self.assertNotIn("as documented", src)
            self.assertIn("interpolated", src)

    def test_pinned_mass_accuracy_is_never_re_measured(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "Orbitrap Exploris 480", 120000, 30000)
            s = dp.parallel_safe(cfg)
            self.assertTrue(s["ok"])
            self.assertEqual(s["measure"], ["window"])

    def test_router_and_generators_real_gate_agree(self):
        """As upstream's RouterAndGeneratorAgreeTests: the generator is RUN (a subprocess, as
        run_search runs it), and the expected answer is pinned too."""
        cases = (("planned 120k/15k", ("Orbitrap Exploris 480", 120000, 15000), None, True),
                 ("pinned 120k/30k", ("Orbitrap Exploris 480", 120000, 30000), None, True),
                 # resolution unknown -> PLAN_AUTO, never measured (see MEASURE_CLASSES): both
                 # levels would be measured and a measured MS1 is what DIA-NN warns about
                 ("auto, resolution unknown", ("Orbitrap Fusion Lumos", None, None), None, False),
                 ("planned, window pinned", ("Orbitrap Exploris 480", 120000, 15000),
                  "--window 7\n", True),
                 ("unknown instrument", ("Mystery-9000", None, None), None, False),
                 ("planned, sidecar gone", ("Orbitrap Exploris 480", 120000, 15000), "RM", False),
                 ("planned, window junk", ("Orbitrap Exploris 480", 120000, 15000),
                  "--window wide\n", False),
                 ("planned, one flag", ("Orbitrap Exploris 480", 120000, 15000),
                  "--mass-acc-ms1 7\n", False))
        for label, (instr, r1, r2), extra, expected in cases:
            with self.subTest(label), tempfile.TemporaryDirectory() as d:
                cfg, _ = estimate(d, instr, r1, r2)
                if extra == "RM":
                    os.remove(cfg + ".rationale.json")
                elif extra:
                    with open(cfg, "a") as fh:
                        fh.write(extra)
                routed, why = run_search.parallel_decision("diann", FILES, cfg, Args())
                raws, fasta = _inputs(d)
                p = _generate(d, cfg, raws, fasta)
                self.assertNotIn("Traceback", p.stderr)
                self.assertEqual(routed, p.returncode == 0, f"{why} | {p.stderr.strip()}")
                self.assertEqual(routed, expected, why)

    def test_the_chains_search_provenance_records_mass_accuracy_as_measured_at_run_time(self):
        """run_search.py end to end on a SLURM host (sbatch stubbed): search_provenance.json
        must not describe a mass accuracy step 1b measures as "not in the cfg -- DIA-NN
        calibrates it itself"."""
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "Orbitrap Exploris 480", 120000, 15000)
            bindir = os.path.join(d, "bin")
            os.makedirs(bindir)
            with open(os.path.join(bindir, "sbatch"), "w") as fh:
                fh.write("#!/bin/sh\necho 1\n")
            os.chmod(os.path.join(bindir, "sbatch"), 0o755)
            raws, fasta = _inputs(d)
            tools, bundle = os.path.join(d, "tools.json"), os.path.join(d, "bundle.json")
            with open(tools, "w") as fh:
                json.dump({"diann": "/bin/true"}, fh)
            with open(bundle, "w") as fh:
                json.dump({"acquisition": "DIA"}, fh)
            env = {k: v for k, v in os.environ.items() if k != "SLURM_JOB_ID"}
            env["PATH"] = bindir + os.pathsep + env.get("PATH", "")
            out = os.path.join(d, "out")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                                "--tools", tools, "--bundle", bundle, "--params", cfg,
                                "--fasta", fasta, "--out", out, "--files", *raws],
                               cwd=d, env=env, capture_output=True, text=True, timeout=120)
            self.assertEqual(p.returncode, 0, p.stderr)
            prov = json.load(open(os.path.join(out, "search_provenance.json")))
            self.assertEqual(prov["search_mode"], "parallel_5step")
            ma = prov["result"]["mass_acc"]
            self.assertTrue(ma["measured"], ma)
            self.assertEqual(ma["value_file"], os.path.join(out, "massacc.txt"))
            self.assertNotIn("calibrates it itself", json.dumps(ma))
            self.assertEqual(prov["result"]["step1b_measures"], ["window", "mass-acc"])
            self.assertEqual(prov["resolved_params_produced"], "runtime")
            self.assertIn("mass accuracy", prov["resolved_params_note"])


class ChangelogTests(unittest.TestCase):
    """A one-flag mass-accuracy override on an instrument with no table value for the other level
    used to produce a cfg and now exits non-zero. That is the right call -- written alone the flag
    fixes the other level at 20 ppm -- but it breaks `workflow.yaml` files that are in use, and a
    behaviour change a user cannot see coming belongs in the changelog."""

    PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(HERE))), "CHANGELOG.md")

    def test_the_breaking_override_change_is_in_the_changelog(self):
        text = open(self.PATH).read()
        self.assertIn("BREAKING", text)
        head = text[:text.index("\n## ", 20)]          # the newest entry only
        self.assertIn("param_overrides", head)
        self.assertIn("--mass-acc", head)
        self.assertIn("20 ppm", head, "the changelog does not say what went wrong before")
        for remedy in ("both", "--ms1-resolution"):
            self.assertIn(remedy, head, f"no remedy: {remedy}")

if __name__ == "__main__":
    unittest.main(verbosity=2)
