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


if __name__ == "__main__":
    unittest.main(verbosity=2)
