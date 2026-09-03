#!/usr/bin/env python3
"""
An unpinned --window must not silently demote a big cohort to a sequential search.

`estimate_params.py` deliberately omits --window: the scan-window radius depends on the
acquisition scheme (cycle time vs peak width) and has to be MEASURED, not guessed. The
5-step chain handles that itself -- step 1b runs probe_window.py and pins one radius into
steps 2 and 4 -- and `diann_parallel.mass_acc_status` is the shared definition of
"parallel-safe", which reports such a cfg as NOT fixed.

So the router must apply the same recoverability rule the chain applies. Reading
mass_acc_status literally made run_search decline and fall back to ONE single-shot search.
That is not a slow path, it is a different order of magnitude: --threads parallelises within
a run, not across runs, so at ~30 min/file a 310-file cohort is ~155 h sequential against a
few hours for the chain. Nothing errors; the user just waits a week.

Observed on a real 310-file timsTOF cohort (PROT_0793, 2026-09-03) with mass accuracy
correctly pinned at 15/15:

    [run_search] parallel routing: no -- 310 files, but mass accuracy is not pinned in
    wf/params_mouse.cfg (not set: --window)
"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import run_search  # noqa: E402
import diann_parallel as dp  # noqa: E402


class Args:
    """The subset of the argparse namespace parallel_decision reads."""
    def __init__(self, threshold=5, no_parallel=False):
        self.parallel_threshold = threshold
        self.no_parallel = no_parallel


def _cfg(tmpdir, body):
    p = os.path.join(tmpdir, "params.cfg")
    with open(p, "w") as fh:
        fh.write(body)
    return p


class WindowIsRecoverableTests(unittest.TestCase):
    """SLURM presence is stubbed, not skipped on. The routing rule under test is about the
    cfg, not about the host: skipping wherever `sbatch` is absent would silence this test on
    every developer laptop and in CI -- i.e. everywhere except the one machine that already
    works."""

    def setUp(self):
        self._real = run_search.slurm_available
        run_search.slurm_available = lambda: True

    def tearDown(self):
        run_search.slurm_available = self._real

    def test_pinned_mass_accuracy_without_window_still_routes_parallel(self):
        """The exact cfg estimate_params.py emits for a timsTOF run."""
        with tempfile.TemporaryDirectory() as d:
            cfg = _cfg(d, "--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n")
            use, why = run_search.parallel_decision("diann", ["f%d.d" % i for i in range(310)],
                                                    cfg, Args())
            self.assertTrue(use, f"declined a 310-file chain over an unpinned --window: {why}")
            self.assertIn("window", why)

    def test_unpinned_mass_accuracy_still_declines(self):
        """Mass accuracy is NOT recoverable -- DIA-NN calibrates it per run against the
        library, so there is no single value step 1b could carry forward."""
        with tempfile.TemporaryDirectory() as d:
            cfg = _cfg(d, "--qvalue 0.01\n")
            use, why = run_search.parallel_decision("diann", ["f%d.d" % i for i in range(310)],
                                                    cfg, Args())
            self.assertFalse(use)
            self.assertIn("mass accuracy", why)

    def test_auto_zero_mass_accuracy_still_declines(self):
        """0 means 'optimise automatically', which is as unsafe as omitting it."""
        with tempfile.TemporaryDirectory() as d:
            cfg = _cfg(d, "--mass-acc 0\n--mass-acc-ms1 0\n")
            use, why = run_search.parallel_decision("diann", ["f%d.d" % i for i in range(310)],
                                                    cfg, Args())
            self.assertFalse(use)

    def test_fully_pinned_cfg_routes_parallel(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = _cfg(d, "--mass-acc 15\n--mass-acc-ms1 15\n--window 7\n")
            use, why = run_search.parallel_decision("diann", ["f%d.d" % i for i in range(310)],
                                                    cfg, Args())
            self.assertTrue(use)


class WindowZeroAndJunkTests(unittest.TestCase):
    """`--window 0` is recoverable; `--window <junk>` is not; and in neither case may a
    --window from the cfg reach the same command line as the measured one."""

    def setUp(self):
        self._real = run_search.slurm_available
        run_search.slurm_available = lambda: True

    def tearDown(self):
        run_search.slurm_available = self._real

    def test_window_zero_routes_and_cfg_window_is_dropped(self):
        """0 is DIA-NN for "optimise automatically" -- per file, which is exactly the
        inconsistency step 1b removes, so it is a case to FIX, not a case to decline.

        But routing it is only safe if the cfg's own --window comes out: steps 2 and 4 get
        `--window $(cat window.txt)` PREFIXED, so a surviving `--window 0` puts two values on
        one line. DIA-NN then either honours the wrong one (silently undoing step 1b, no
        error) or rejects 0 outright and kills every array task, after step 1 has already
        paid for the library prediction.
        """
        with tempfile.TemporaryDirectory() as d:
            cfg = _cfg(d, "--mass-acc 15\n--mass-acc-ms1 15\n--window 0\n--qvalue 0.01\n")
            use, why = run_search.parallel_decision("diann", ["f%d.d" % i for i in range(310)],
                                                    cfg, Args())
            self.assertTrue(use, f"declined a recoverable --window 0: {why}")
            flags = dp.read_cfg_flags(cfg, drop=("--window",))
            self.assertNotIn("--window", flags,
                             "cfg --window survived into the step flags; it would collide "
                             "with the measured one")
            self.assertIn("--qvalue", flags, "dropping --window must not drop anything else")

    def test_unparseable_window_declines(self):
        """A typo is a mistake to report, not something to quietly measure over."""
        with tempfile.TemporaryDirectory() as d:
            cfg = _cfg(d, "--mass-acc 15\n--mass-acc-ms1 15\n--window wide\n")
            use, why = run_search.parallel_decision("diann", ["f%d.d" % i for i in range(310)],
                                                    cfg, Args())
            self.assertFalse(use, "measured over a typo'd --window instead of reporting it")
            self.assertIn("window", why)

    def test_negative_mass_accuracy_declines(self):
        """`not in (None, 0)` accepted any non-zero float, so a corrupt cfg routed a whole
        cohort to the cluster; DIA-NN exits 0 on fatal errors, so it surfaces late."""
        with tempfile.TemporaryDirectory() as d:
            cfg = _cfg(d, "--mass-acc -3\n--mass-acc-ms1 15\n")
            use, why = run_search.parallel_decision("diann", ["f%d.d" % i for i in range(310)],
                                                    cfg, Args())
            self.assertFalse(use)


class RouterAndGeneratorAgreeTests(unittest.TestCase):
    """The original bug was DRIFT: run_search decided one way, diann_parallel the other, and
    the router won. Pinning each side's answer separately would not have caught it -- this
    asserts they are the SAME answer, for every shape of cfg we know about."""

    CASES = [
        ("pinned, window unset", "--mass-acc 15\n--mass-acc-ms1 15\n"),
        ("pinned, window 0", "--mass-acc 15\n--mass-acc-ms1 15\n--window 0\n"),
        ("fully pinned", "--mass-acc 15\n--mass-acc-ms1 15\n--window 7\n"),
        ("no mass accuracy", "--qvalue 0.01\n"),
        ("mass accuracy auto", "--mass-acc 0\n--mass-acc-ms1 0\n"),
        ("mass accuracy negative", "--mass-acc -3\n--mass-acc-ms1 15\n"),
        ("window junk", "--mass-acc 15\n--mass-acc-ms1 15\n--window wide\n"),
        ("ms1 only", "--mass-acc-ms1 15\n"),
    ]

    def setUp(self):
        self._real = run_search.slurm_available
        run_search.slurm_available = lambda: True

    def tearDown(self):
        run_search.slurm_available = self._real

    def test_router_matches_generator_for_every_cfg_shape(self):
        with tempfile.TemporaryDirectory() as d:
            for label, body in self.CASES:
                cfg = os.path.join(d, label.replace(" ", "_").replace(",", "") + ".cfg")
                with open(cfg, "w") as fh:
                    fh.write(body)
                routed, why = run_search.parallel_decision(
                    "diann", ["f%d.d" % i for i in range(310)], cfg, Args())
                # the generator's own verdict, with the defaults run_diann_parallel passes
                generated = dp.parallel_safe(cfg)["ok"]
                self.assertEqual(routed, generated,
                                 f"{label}: router says {routed}, generator says "
                                 f"{generated} -- they must not disagree ({why})")


if __name__ == "__main__":
    unittest.main(verbosity=2)
