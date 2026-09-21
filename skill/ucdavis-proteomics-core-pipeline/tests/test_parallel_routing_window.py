#!/usr/bin/env python3
"""
An unpinned --window must not silently demote a big cohort to a sequential search -- and the
chain it routes to must actually run what the gate approved.

`estimate_params.py` deliberately omits --window: the scan-window radius depends on the
acquisition scheme (cycle time vs peak width) and has to be MEASURED, not guessed. The
5-step chain handles that itself -- step 1b runs probe_window.py and pins one radius into
steps 2-5 -- and `diann_parallel.parallel_safe` is the one definition of "may this cfg run
as the chain", consumed by both the router (run_search.py) and the generator.

Reading the window as "not parallel-safe" made run_search decline and fall back to ONE
single-shot search. That is not a slow path, it is a different order of magnitude: --threads
parallelises within a run, not across runs, so at ~30 min/file a 310-file cohort is ~155 h
sequential against a few hours for the chain. Nothing errors; the user just waits a week.

Observed on a real 310-file timsTOF cohort (PROT_0793, 2026-09-03) with mass accuracy
correctly pinned at 15/15:

    [run_search] parallel routing: no -- 310 files, but mass accuracy is not pinned in
    wf/params_mouse.cfg (not set: --window)

Routing those cohorts to the chain made step 1b the common path, so the review of that fix
is covered here too: step 1b must get the .NET 8 environment for Thermo .raw, survive a blank
first file, and never leave a "resolved" cfg without a window; the gate and the generated
steps must read the cfg by one rule; and `--sbatch job.sh && sbatch job.sh` must not be able
to resubmit a stale job script.
"""
import contextlib
import glob
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import run_search  # noqa: E402
import diann_parallel as dp  # noqa: E402
import estimate_params  # noqa: E402

GENERATOR = os.path.join(SCRIPTS, "diann_parallel.py")
PINNED_MA = "--mass-acc 15\n--mass-acc-ms1 15\n"


class Args:
    """The subset of the argparse namespace parallel_decision reads."""
    def __init__(self, threshold=5, no_parallel=False):
        self.parallel_threshold = threshold
        self.no_parallel = no_parallel


def _cfg(tmpdir, body, name="params.cfg"):
    p = os.path.join(tmpdir, name)
    with open(p, "w") as fh:
        fh.write(body)
    return p


def _read(path):
    with open(path) as fh:
        return fh.read()


def _write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def _route(cfg, n=310):
    return run_search.parallel_decision("diann", ["f%d.d" % i for i in range(n)], cfg, Args())


def _inputs(d, names=("f0.d", "f1.d", "f2.d", "f3.d", "f4.d", "f5.d")):
    """Real paths for a cohort: .d are directories, anything else a file."""
    raws = []
    for name in names:
        p = os.path.join(d, name)
        if name.endswith(".d"):
            os.makedirs(p, exist_ok=True)
        else:
            open(p, "w").close()
        raws.append(p)
    fasta = os.path.join(d, "db.fasta")
    with open(fasta, "w") as fh:
        fh.write(">sp|P1|X\nPEPTIDER\n")
    return raws, fasta


def _generate(d, cfg, raws, fasta, *extra):
    """Run the generator exactly as run_search.run_diann_parallel does: a subprocess."""
    out = os.path.join(d, "out")
    p = subprocess.run([sys.executable, GENERATOR, "--diann", "/bin/true", "--raw", *raws,
                        "--fasta", fasta, "--out", out, "--cfg", cfg, *extra],
                       capture_output=True, text=True)
    return p, out


class _SlurmStub(unittest.TestCase):
    """SLURM presence is stubbed, not skipped on. The routing rule under test is about the
    cfg, not about the host: skipping wherever `sbatch` is absent would silence these tests on
    every developer laptop and in CI -- i.e. everywhere except the one machine that already
    works."""

    def setUp(self):
        self._real = run_search.slurm_available
        run_search.slurm_available = lambda: True

    def tearDown(self):
        run_search.slurm_available = self._real


class WindowIsRecoverableTests(_SlurmStub):

    def test_pinned_mass_accuracy_without_window_still_routes_parallel(self):
        """The exact cfg estimate_params.py emits for a timsTOF run."""
        with tempfile.TemporaryDirectory() as d:
            use, why = _route(_cfg(d, "--qvalue 0.01\n" + PINNED_MA))
            self.assertTrue(use, f"declined a 310-file chain over an unpinned --window: {why}")
            self.assertIn("window", why)

    def test_unpinned_mass_accuracy_still_declines(self):
        """Unpinned mass accuracy with no estimate_params.py plan to measure it declines. (The
        one recoverable case -- an Orbitrap with no documented tier, measured in step 1b -- is
        in test_orbitrap_mass_accuracy.py; this cfg has no sidecar, so it is not that case.)"""
        with tempfile.TemporaryDirectory() as d:
            use, why = _route(_cfg(d, "--qvalue 0.01\n"))
            self.assertFalse(use)
            self.assertIn("mass accuracy", why)

    def test_zero_mass_accuracy_declines_and_is_not_called_auto(self):
        """0 is NOT auto in DIA-NN: it is a literal 0 ppm tolerance ("Mass accuracy will be
        fixed to 0 (MS2) and 0 (MS1)") and a 28-run Lumos search returned 0 IDs (PR #38).
        The reason used to say "auto (0)", which sends the reader the wrong way."""
        with tempfile.TemporaryDirectory() as d:
            use, why = _route(_cfg(d, "--mass-acc 0\n--mass-acc-ms1 0\n"))
            self.assertFalse(use)
            self.assertIn("0 ppm", why)
            self.assertNotIn("auto (0)", why)

    def test_fully_pinned_cfg_routes_parallel(self):
        with tempfile.TemporaryDirectory() as d:
            use, _ = _route(_cfg(d, PINNED_MA + "--window 7\n"))
            self.assertTrue(use)


class WindowZeroAndJunkTests(_SlurmStub):
    """`--window 0` is recoverable; anything that is not a non-negative integer is not; and in
    neither case may a --window from the cfg reach the same command line as the measured one."""

    def test_window_zero_routes_and_cfg_window_is_dropped(self):
        """DIA-NN does not accept 0 as a radius ("scan window radius should be a positive
        integer") and optimises per file instead -- seen on the 18-file poplar run -- which is
        exactly the inconsistency step 1b removes, so it is a case to FIX, not to decline.

        But routing it is only safe if the cfg's own --window comes out: steps 2-5 get
        `--window $(cat window.txt)` PREFIXED, so a surviving `--window 0` puts two values on
        one line."""
        with tempfile.TemporaryDirectory() as d:
            cfg = _cfg(d, PINNED_MA + "--window 0\n--qvalue 0.01\n")
            use, why = _route(cfg)
            self.assertTrue(use, f"declined a recoverable --window 0: {why}")
            flags = dp.read_cfg_flags(cfg, drop=("--window",))
            self.assertNotIn("--window", flags,
                             "cfg --window survived into the step flags; it would collide "
                             "with the measured one")
            self.assertIn("--qvalue", flags, "dropping --window must not drop anything else")

    def test_unparseable_window_declines(self):
        """A typo is a mistake to report, not something to quietly measure over."""
        with tempfile.TemporaryDirectory() as d:
            use, why = _route(_cfg(d, PINNED_MA + "--window wide\n"))
            self.assertFalse(use, "measured over a typo'd --window instead of reporting it")
            self.assertIn("window", why)

    def test_negative_mass_accuracy_declines(self):
        """`not in (None, 0)` accepted any non-zero float, so a corrupt cfg routed a whole
        cohort to the cluster; DIA-NN exits 0 on fatal errors, so it surfaces late."""
        with tempfile.TemporaryDirectory() as d:
            use, _ = _route(_cfg(d, "--mass-acc -3\n--mass-acc-ms1 15\n"))
            self.assertFalse(use)

    def test_window_must_be_a_positive_integer(self):
        """Finding 10. DIA-NN wants a positive INTEGER. `--window 0.5` read as a truthy float,
        skipped the probe and went straight to DIA-NN; `--window nan` crashed int() with a
        traceback inside the router."""
        for bad in ("0.5", "nan", "inf", "-1", "7.0", "wide", "1e3"):
            with self.subTest(window=bad), tempfile.TemporaryDirectory() as d:
                cfg = _cfg(d, PINNED_MA + f"--window {bad}\n")
                safe = dp.parallel_safe(cfg)
                self.assertFalse(safe["ok"], f"--window {bad} was accepted")
                self.assertEqual(safe["code"], "window_invalid")
                self.assertIn("positive integer", safe["reason"])
                use, why = _route(cfg)          # must not raise
                self.assertFalse(use)
        with tempfile.TemporaryDirectory() as d:
            safe = dp.parallel_safe(_cfg(d, PINNED_MA + "--window 7\n"))
            self.assertEqual((safe["ok"], safe["probe"], safe["code"]), (True, False, "pinned"))


class OneCfgRuleTests(_SlurmStub):
    """Findings 2, 3, 11. The gate read the cfg by shlex TOKEN, the step flags by LINE, and
    params.base.cfg by `split()[:1]`. Where they disagreed, the gate approved flags the
    generated steps did not actually carry."""

    def test_a_midline_or_tab_separated_window_is_dropped_everywhere(self):
        for label, body in (("mid-line", "--mass-acc 15 --mass-acc-ms1 15 --window 0 --qvalue 0.01\n"),
                            ("tab", PINNED_MA + "--window\t0\n--qvalue 0.01\n")):
            with self.subTest(label), tempfile.TemporaryDirectory() as d:
                cfg = _cfg(d, body)
                safe = dp.parallel_safe(cfg)
                self.assertTrue(safe["ok"] and safe["probe"], safe["reason"])
                self.assertNotIn("--window", dp.read_cfg_flags(cfg, drop=("--window",)))
                raws, fasta = _inputs(d)
                p, out = _generate(d, cfg, raws, fasta)
                self.assertEqual(p.returncode, 0, p.stderr)
                for step in ("step2_firstpass.sbatch", "step3_assembly.sbatch",
                             "step4_finalpass.sbatch", "step5_report.sbatch"):
                    body_s = _read(os.path.join(out, step))
                    self.assertEqual(body_s.count("--window"), 1,
                                     f"{label}: {step} carries more than the measured --window")
                    self.assertIn("--window $(cat", body_s)
                    self.assertIn("--qvalue 0.01", body_s)
                base = _read(os.path.join(out, "params.base.cfg"))
                self.assertNotIn("--window", base, f"{label}: params.base.cfg kept the cfg --window")
                self.assertIn("--qvalue 0.01", base)

    def test_an_inline_comment_does_not_swallow_the_flags_after_it(self):
        """`--qvalue 0.01   # 1% FDR` was spliced into one joined bash line, so bash treated
        every later flag as part of the comment -- while the gate, reading tokens, counted
        them as set."""
        with tempfile.TemporaryDirectory() as d:
            cfg = _cfg(d, "--qvalue 0.01   # 1% FDR\n--mass-acc 15\n--mass-acc-ms1 15\n"
                          "--window 7   # measured on PROT_0793\n--cut K*,R*\n")
            self.assertTrue(dp.parallel_safe(cfg)["ok"])
            flags = dp.read_cfg_flags(cfg)
            argv = subprocess.run(["bash", "-c", f"printf '%s\\n' {flags}"],
                                  capture_output=True, text=True, cwd=d).stdout.split("\n")
            self.assertNotIn("#", flags)
            for want in ("--mass-acc", "15", "--mass-acc-ms1", "--window", "7", "K*,R*"):
                self.assertIn(want, argv, f"{want} did not reach DIA-NN's command line")

    def test_a_quoted_value_with_spaces_stays_one_word(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = _cfg(d, PINNED_MA + '--lib-dir "/data/my libs"\n')
            flags = dp.read_cfg_flags(cfg)
            argv = subprocess.run(["bash", "-c", f"printf '%s\\n' {flags}"],
                                  capture_output=True, text=True).stdout.split("\n")
            self.assertIn("/data/my libs", argv)
            # and the base cfg reads back to the same tokens by the same rule
            base = os.path.join(d, "base.cfg")
            dp.write_cfg(cfg, base, drop=("--window",))
            self.assertEqual(dp.cfg_tokens(base), dp.cfg_tokens(cfg))


class RouterAndGeneratorAgreeTests(_SlurmStub):
    """Findings 4 and 5. The original bug was DRIFT: run_search decided one way, diann_parallel
    the other. The earlier version of this test compared parallel_safe() with itself, so it
    could not see the generator's REAL gate, which read ma["fixed"] -- and a duplicated,
    unparseable --mass-acc passed that gate. This runs the generator the way run_search does
    (a subprocess) and pins the expected answer too, so both sides flipping together is caught.
    """

    CASES = [
        ("pinned, window unset", PINNED_MA, True),
        ("pinned, window 0", PINNED_MA + "--window 0\n", True),
        ("pinned, window 0 mid-line", "--mass-acc 15 --mass-acc-ms1 15 --window 0\n", True),
        ("fully pinned", PINNED_MA + "--window 7\n", True),
        ("no mass accuracy", "--qvalue 0.01\n", False),
        ("mass accuracy 0", "--mass-acc 0\n--mass-acc-ms1 0\n", False),
        ("mass accuracy negative", "--mass-acc -3\n--mass-acc-ms1 15\n", False),
        ("mass accuracy inf", "--mass-acc inf\n--mass-acc-ms1 15\n", False),
        ("duplicate bad mass accuracy", PINNED_MA + "--mass-acc abc\n", False),
        ("conflicting mass accuracy", PINNED_MA + "--mass-acc 10\n", False),
        ("window junk", PINNED_MA + "--window wide\n", False),
        ("window fractional", PINNED_MA + "--window 0.5\n", False),
        ("window nan", PINNED_MA + "--window nan\n", False),
        ("ms1 only", "--mass-acc-ms1 15\n", False),
        ("unbalanced quote", PINNED_MA + '--var-mod "UniMod:35,15.994915,M\n', False),
    ]

    def test_router_matches_the_generators_real_gate_for_every_cfg_shape(self):
        for label, body, expected in self.CASES:
            with self.subTest(label), tempfile.TemporaryDirectory() as d:
                cfg = _cfg(d, body)
                routed, why = _route(cfg)
                raws, fasta = _inputs(d)
                p, _ = _generate(d, cfg, raws, fasta)
                generated = p.returncode == 0
                self.assertNotIn("Traceback", p.stderr, f"{label}: generator crashed")
                self.assertEqual(routed, generated,
                                 f"{label}: router says {routed}, generator says {generated} "
                                 f"-- they must not disagree ({why} | {p.stderr.strip()})")
                self.assertEqual(routed, expected, f"{label}: {why}")

    def test_the_override_is_only_for_omitted_mass_accuracy(self):
        """--allow-auto-mass-acc is named for, and limited to, mass accuracy left to DIA-NN.
        It must not wave through a duplicated junk value, which parallel_safe used to hide on
        the probe path by overwriting fixed=True."""
        with tempfile.TemporaryDirectory() as d:
            raws, fasta = _inputs(d)
            p, _ = _generate(d, _cfg(d, PINNED_MA + "--mass-acc abc\n"), raws, fasta,
                             "--allow-auto-mass-acc")
            self.assertNotEqual(p.returncode, 0, "an invalid mass accuracy was overridden")
            self.assertIn("not a number", p.stderr)
            # omitted mass accuracy must not carry a junk --window through the override
            p, _ = _generate(d, _cfg(d, "--window wide\n", "junkwin.cfg"), raws, fasta,
                             "--allow-auto-mass-acc")
            self.assertNotEqual(p.returncode, 0, "a junk --window rode through the override")
            p, _ = _generate(d, _cfg(d, "--qvalue 0.01\n", "auto.cfg"), raws, fasta,
                             "--allow-auto-mass-acc")
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("WARNING", p.stderr)


class RemedyFollowsTheReasonTests(_SlurmStub):
    """Finding 9. Every decline said "re-run estimate_params.py with the instrument table" --
    including a typo'd --window, which estimate_params.py never writes -- and the generator told
    users to hand-run probe_window.py for a window step 1b measures by itself."""

    def test_the_fix_named_matches_the_cause(self):
        table = estimate_params.instrument_ppm_summary()
        with tempfile.TemporaryDirectory() as d:
            _, why = _route(_cfg(d, PINNED_MA + "--window wide\n", "w.cfg"))
            self.assertNotIn("estimate_params", why)
            self.assertIn("step 1b", why)
            _, why = _route(_cfg(d, "--qvalue 0.01\n", "m.cfg"))
            self.assertIn("estimate_params.py", why)
            self.assertIn(table, why, "the ppm table text is not the one estimate_params owns")

    def test_the_generator_does_not_send_users_to_probe_by_hand(self):
        with tempfile.TemporaryDirectory() as d:
            raws, fasta = _inputs(d)
            for body in ("--qvalue 0.01\n", PINNED_MA + "--window wide\n"):
                p, _ = _generate(d, _cfg(d, body), raws, fasta)
                self.assertNotEqual(p.returncode, 0)
                self.assertNotIn("probe_window", p.stderr)

    def test_the_table_is_read_from_classify_instruments_own_rows(self):
        for name, key in (("timsTOF HT", "timstof"), ("Orbitrap Astral", "orbitrap_astral"),
                          ("ZenoTOF 7600", "sciex_tof")):
            _, ms1, ms2, _, _ = estimate_params.classify_instrument(name)
            self.assertEqual((ms1, ms2), estimate_params.DIANN_INSTRUMENT_PPM[key])
            self.assertIn(f"{ms1}/{ms2}", estimate_params.instrument_ppm_summary())


# Stand-in for DIA-NN 2.6. Like the real binary it exits 0 on a fatal error; like the real
# binary it cannot read a Thermo .raw without .NET 8 in its environment.
FAKE_DIANN = r"""#!/bin/bash
raw=""
while [ $# -gt 0 ]; do
  if [ "$1" = "--f" ]; then raw="$2"; fi
  shift
done
echo "$(basename "$raw")" >> "$FAKE_DIANN_LOG"
case "$raw" in
  *.raw) if [ -z "$DOTNET_ROOT" ]; then echo "ERROR: cannot open $raw: .NET runtime not found"; exit 0; fi ;;
esac
case "$(basename "$raw")" in
  blank*) echo "WARNING: no precursors identified, calibration skipped"; exit 0 ;;
esac
echo "Scan window radius set to 7"
"""

FAKE_DOTNET = 'export DOTNET_ROOT=/opt/fake-dotnet8; export PATH=/opt/fake-dotnet8:"$PATH"; '


class Step1bRunsTests(unittest.TestCase):
    """Findings 1, 6, 12 -- by EXECUTING the generated step1b_window.sbatch against a fake
    DIA-NN, not by reading its text."""

    def _chain(self, d, names):
        raws, fasta = _inputs(d, names)
        diann = os.path.join(d, "fake-diann")
        _write(diann, FAKE_DIANN)
        os.chmod(diann, 0o755)
        cfg = _cfg(d, "--qvalue 0.01\n" + PINNED_MA)
        out = os.path.join(d, "out")
        argv = [GENERATOR, "--diann", diann, "--raw", *raws, "--fasta", fasta, "--out", out,
                "--cfg", cfg]
        # ensure_dotnet8.sh installs a real runtime over the network; stub only that. The .raw
        # detection and everything the prefix is spliced into are the real code.
        fake_prefix = lambda rs: FAKE_DOTNET if any(r.endswith(".raw") for r in rs) else ""
        with mock.patch.object(sys, "argv", argv), \
                mock.patch.object(dp, "dotnet_prefix", fake_prefix), \
                contextlib.redirect_stdout(io.StringIO()) as so:
            dp.main()
        info = json.loads(so.getvalue())
        _write(os.path.join(out, "step1.predicted.speclib"), "lib")   # step 1's artefact
        return out, info

    def _run_step1b(self, d, out):
        log = os.path.join(d, "probed.txt")
        env = {k: v for k, v in os.environ.items() if k != "DOTNET_ROOT"}
        env["FAKE_DIANN_LOG"] = log
        p = subprocess.run(["bash", os.path.join(out, "step1b_window.sbatch")], cwd=out,
                           env=env, capture_output=True, text=True, timeout=120)
        probed = _read(log).split() if os.path.exists(log) else []
        return p, probed

    def test_thermo_raw_probe_gets_the_dotnet_environment(self):
        """Finding 1. Step 1b passed plain --diann while every other step ran the .NET-8
        prefixed command, so on a .raw cohort DIA-NN could not read the file, window.txt stayed
        empty and steps 2-5 waited on afterok for ever."""
        with tempfile.TemporaryDirectory() as d:
            out, _ = self._chain(d, ["s%d.raw" % i for i in range(6)])
            p, probed = self._run_step1b(d, out)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(_read(os.path.join(out, "window.txt")).strip(), "7")
            # every representative run answered, so none was replaced
            self.assertEqual(len(probed), dp.PROBE_CANDIDATES, probed)
            self.assertTrue(all(r.endswith(".raw") for r in probed), probed)

    def test_a_run_that_gives_no_radius_falls_through_to_the_next_representative_run(self):
        """Finding 12. Step 1b probed raws[0] only, so one blank or wash first in the list
        failed the whole cohort. The radius is a property of the method, so a run that gives no
        radius is replaced -- by the next run nearest the median (tests/test_step1b_window_probe.py
        covers the selection), not by the next file of the listing. The files here are all
        empty, so the blank cannot be told apart by size and lands on the median itself."""
        with tempfile.TemporaryDirectory() as d:
            out, info = self._chain(d, ["a0.mzML", "a1.mzML", "blank_01.mzML", "c0.mzML",
                                        "c1.mzML"])
            self.assertEqual(info["scan_window"]["evidence_file"],
                             os.path.join(out, "window.json"))
            p, probed = self._run_step1b(d, out)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(probed, ["blank_01.mzML", "a1.mzML", "c0.mzML", "c1.mzML"])
            w = json.loads(_read(os.path.join(out, "window.json")))
            self.assertEqual(w["failed"], ["blank_01.mzML"])
            self.assertEqual(w["radii"], [7, 7, 7])
            resolved = dp.cfg_tokens(os.path.join(out, "params.resolved.cfg"))
            self.assertEqual(resolved.count("--window"), 1)
            self.assertEqual(resolved[resolved.index("--window") + 1], "7")

    def test_no_radius_fails_loudly_and_leaves_no_resolved_cfg(self):
        """Finding 6. params.resolved.cfg was copied into place BEFORE the probe ran, so a failed
        probe left a "resolved" cfg with no window -- and a resubmit could find the previous
        run's window.txt. Both must be gone, and the job must fail."""
        with tempfile.TemporaryDirectory() as d:
            out, info = self._chain(d, ["blank_%d.mzML" % i for i in range(6)])
            self.assertEqual(info["resolved_params"]["produced"], "runtime")
            self.assertFalse(os.path.exists(os.path.join(out, "params.resolved.cfg")),
                             "the generator created the resolved cfg before anything was measured")
            for stale in ("window.txt", "params.resolved.cfg"):           # a previous run's
                _write(os.path.join(out, stale), "--window 9\n")
            p, probed = self._run_step1b(d, out)
            self.assertNotEqual(p.returncode, 0)
            self.assertEqual(len(probed), dp.PROBE_MAX_FAILURES)
            self.assertIn("FAILED", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(out, "window.txt")))
            self.assertFalse(os.path.exists(os.path.join(out, "params.resolved.cfg")))
            # Round 2, item 5: resubmitting step 1b alone does not restart the chain -- steps
            # 2-5 are afterok on THIS job id. The message must say so and name the way out.
            self.assertIn("DependencyNeverSatisfied", p.stderr)
            self.assertIn("jobs.txt", p.stderr)
            self.assertIn("steps 2-5", p.stderr)

    def test_a_resolved_cfg_that_cannot_be_written_says_so_not_dianns_fault(self):
        """Round 2, item 5. The check on params.resolved.cfg reused must_exist(), whose message
        is "DIA-NN exited 0 but did not write ..." -- but bash writes that file, not DIA-NN.
        A directory of that name makes the move land INSIDE it, which `-s` alone accepts."""
        with tempfile.TemporaryDirectory() as d:
            out, _ = self._chain(d, ["s%d.mzML" % i for i in range(6)])
            blocker = os.path.join(out, "params.resolved.cfg")
            os.makedirs(blocker)
            _write(os.path.join(blocker, "keep"), "x")
            p, _ = self._run_step1b(d, out)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("could not be moved into place", p.stderr)
            self.assertNotIn("DIA-NN exited 0", p.stderr)

    def test_step1b_wall_clock_covers_every_attempt(self):
        """Round 2, item 7. The per-attempt timeout went from probe_window's 3600 s to 2700 s to
        fit three attempts into 3 h -- unmeasured on a large Astral .raw. Keep 3600 s and size
        the wall clock to the attempts instead. The probe's --budget is that wall clock less a
        margin, so replacements after a failure can never run the job into SLURM's limit before
        window.json is written."""
        self.assertEqual(dp.PROBE_TIMEOUT_S, 3600)
        self.assertGreater(dp.PROBE_WALL_HOURS * 3600, dp.PROBE_CANDIDATES * dp.PROBE_TIMEOUT_S)
        self.assertGreaterEqual(dp.PROBE_BUDGET_S, dp.PROBE_CANDIDATES * dp.PROBE_TIMEOUT_S)
        self.assertLess(dp.PROBE_BUDGET_S, dp.PROBE_WALL_HOURS * 3600)
        with tempfile.TemporaryDirectory() as d:
            out, _ = self._chain(d, ["s%d.mzML" % i for i in range(6)])
            body = _read(os.path.join(out, "step1b_window.sbatch"))
            self.assertIn("--timeout 3600", body)
            self.assertIn(f"--budget {dp.PROBE_BUDGET_S}", body)
            self.assertIn(f"#SBATCH --time={dp.PROBE_WALL_HOURS}:00:00", body)


class ProbeTimeoutTests(unittest.TestCase):
    """probe_window's deadline was only checked when DIA-NN printed a line. A DIA-NN that went
    silent blocked the read until SLURM killed step 1b -- so the fallback to the next file
    (finding 12) would never get its turn."""

    def test_a_silent_diann_is_cut_off_at_the_timeout(self):
        import probe_window
        import time
        with tempfile.TemporaryDirectory() as d:
            fake = os.path.join(d, "silent-diann")
            _write(fake, "#!/bin/bash\nexec sleep 60\n")
            os.chmod(fake, 0o755)
            t0 = time.time()
            radius, _ = probe_window.probe(fake, "x.mzML", "db.fasta", "lib", 1, timeout=1)
            self.assertIsNone(radius)
            self.assertLess(time.time() - t0, 20, "probe hung on a silent DIA-NN")

    def _forking(self, d, body):
        """A --diann that FORKS (no exec), as a bash wrapper does and apptainer very likely
        does: the grandchild holds the stdout pipe. Its pid is written before anything else."""
        pidfile = os.path.join(d, "grandchild.pid")
        fake = os.path.join(d, "forking-diann")
        _write(fake, f"#!/bin/bash\nsleep 300 &\necho $! > {pidfile}\n{body}\nwait\n")
        os.chmod(fake, 0o755)
        return fake, pidfile

    def _gone(self, pidfile):
        """Dead, or a zombie waiting for init to reap it (not burning anything)."""
        import time
        pid = int(_read(pidfile))
        for _ in range(100):
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                return True
            st = subprocess.run(["ps", "-o", "stat=", "-p", str(pid)],
                                capture_output=True, text=True).stdout.strip()
            if not st or st.startswith("Z"):
                return True
            time.sleep(0.1)
        os.kill(pid, 9)                     # do not leak it past the test
        return False

    def test_a_forking_silent_diann_is_cut_off_and_leaves_no_orphan(self):
        """Round 2, item 3. The watchdog killed only the direct child; the forked grandchild kept
        the pipe open, so the read blocked until it exited on its own (60 s here) and kept
        burning the job's CPUs."""
        import probe_window
        import time
        with tempfile.TemporaryDirectory() as d:
            fake, pidfile = self._forking(d, "")
            t0 = time.time()
            radius, _ = probe_window.probe(fake, "x.mzML", "db.fasta", "lib", 1, timeout=1)
            self.assertIsNone(radius)
            # well under the grandchild's 300 s; allows _end_group's 30 s grace on a host whose
            # init is slow to reap
            self.assertLess(time.time() - t0, 60, "probe blocked on a forked grandchild")
            self.assertTrue(self._gone(pidfile), "the forked DIA-NN was left running")

    def test_a_forking_diann_that_answers_leaves_no_orphan(self):
        import probe_window
        with tempfile.TemporaryDirectory() as d:
            fake, pidfile = self._forking(d, 'echo "Scan window radius set to 7"')
            radius, _ = probe_window.probe(fake, "x.mzML", "db.fasta", "lib", 1, timeout=60)
            self.assertEqual(radius, 7)
            self.assertTrue(self._gone(pidfile), "the forked DIA-NN kept running after the answer")


class SbatchCannotResubmitAStaleJobTests(unittest.TestCase):
    """Finding 7. With the chain, `--sbatch job.sh` writes nothing; a NOTE and exit 0 let the
    documented `run_search.py ... --sbatch job.sh && sbatch job.sh` resubmit whatever job.sh was
    already there -- e.g. an old sequential 310-file search."""

    def _run(self, d, existing=None, sbatch="job.sh", cfg=None, before=None):
        bindir = os.path.join(d, "bin")
        os.makedirs(bindir)
        _write(os.path.join(bindir, "sbatch"), "#!/bin/sh\necho 1\n")   # slurm_available() true
        os.chmod(os.path.join(bindir, "sbatch"), 0o755)
        raws, fasta = _inputs(d)
        cfg = cfg or _cfg(d, "--qvalue 0.01\n--xic 10\n--mobilograms\n" + PINNED_MA)
        tools, bundle = os.path.join(d, "tools.json"), os.path.join(d, "bundle.json")
        _write(tools, json.dumps({"diann": "/bin/true"}))
        _write(bundle, json.dumps({"acquisition": "DIA"}))
        if existing is not None:
            _write(os.path.join(d, "job.sh"), existing)
        if before:
            before(d)
        env = {k: v for k, v in os.environ.items() if k != "SLURM_JOB_ID"}
        env["PATH"] = bindir + os.pathsep + env.get("PATH", "")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                            "--tools", tools, "--bundle", bundle, "--params", cfg,
                            "--fasta", fasta, "--out", os.path.join(d, "out"),
                            "--files", *raws, "--sbatch", sbatch],
                           cwd=d, env=env, capture_output=True, text=True)
        return p

    def test_a_stale_job_script_is_moved_aside_and_the_exit_is_nonzero(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._run(d, existing="#!/bin/bash\n# OLD sequential 310-file search\n")
            self.assertEqual(p.returncode, run_search.SBATCH_NOT_WRITTEN, p.stderr)
            self.assertFalse(os.path.exists(os.path.join(d, "job.sh")),
                             "a stale job.sh is still there for `&& sbatch job.sh` to submit")
            moved = glob.glob(os.path.join(d, "job.sh.stale-*"))
            self.assertEqual(len(moved), 1, "the user's file must be renamed, never deleted")
            self.assertIn("OLD sequential", _read(moved[0]))
            self.assertIn("submit.sh", p.stderr)
            self.assertTrue(os.path.exists(os.path.join(d, "out", "submit.sh")))
            prov = json.loads(_read(os.path.join(d, "out", "search_provenance.json")))
            self.assertIsNone(prov["submitted_sbatch"])
            self.assertEqual(os.path.basename(prov["sbatch_refused"]["existing_file_moved_to"]),
                             os.path.basename(moved[0]))
            # finding 6: recorded as produced at run time, not as already resolved
            self.assertEqual(prov["resolved_params_produced"], "runtime")
            self.assertIn("step 1b", prov["scan_window"]["source"])
            # round 2, item 4a: mass accuracy IS pinned; the unset window is not its business
            self.assertTrue(prov["result"]["mass_acc"]["fixed"], prov["result"]["mass_acc"])
            self.assertNotIn("window", prov["result"]["mass_acc"]["reason"])

    def test_sbatch_naming_a_directory_is_refused_and_nothing_moves(self):
        """Round 2, item 2. `--sbatch proj` renamed the project folder -- the cfg inside it too
        -- and generation then failed blaming mass accuracy."""
        with tempfile.TemporaryDirectory() as d:
            proj = os.path.join(d, "proj")
            os.makedirs(proj)
            cfg = _cfg(proj, "--qvalue 0.01\n--xic 10\n--mobilograms\n" + PINNED_MA)
            p = self._run(d, sbatch="proj", cfg=cfg)
            self.assertNotIn(p.returncode, (0, run_search.SBATCH_NOT_WRITTEN), p.stderr)
            self.assertIn("a directory", p.stderr + p.stdout)
            self.assertTrue(os.path.isfile(cfg), "the folder holding the cfg was moved")
            self.assertEqual(glob.glob(os.path.join(d, "proj.stale-*")), [])
            self.assertNotIn("mass accuracy", p.stderr)

    def test_a_failed_generation_leaves_the_job_script_where_it_was(self):
        """Round 2, item 2: set aside only AFTER the chain exists."""
        def block_generation(d):                    # diann_parallel cannot write its file list
            os.makedirs(os.path.join(d, "out", "file_list.txt"))
        with tempfile.TemporaryDirectory() as d:
            p = self._run(d, existing="#!/bin/bash\n# the user's script\n", before=block_generation)
            self.assertNotEqual(p.returncode, 0)
            self.assertNotEqual(p.returncode, run_search.SBATCH_NOT_WRITTEN)
            self.assertEqual(_read(os.path.join(d, "job.sh")), "#!/bin/bash\n# the user's script\n")
            self.assertEqual(glob.glob(os.path.join(d, "job.sh.stale-*")), [])

    def test_a_missing_cfg_is_reported_as_missing(self):
        """Round 2, item 2. cfg_tokens returned [] for a path that was not there, which every
        reader then took for "mass accuracy is not pinned"."""
        with tempfile.TemporaryDirectory() as d:
            gone = os.path.join(d, "moved-away", "params.cfg")
            safe = dp.parallel_safe(gone)
            self.assertEqual(safe["code"], "cfg_missing")
            self.assertIn("cfg not found", safe["reason"])
            with mock.patch.object(run_search, "slurm_available", lambda: True):
                use, why = _route(gone)
            self.assertFalse(use)
            self.assertIn("cfg not found", why)
            self.assertNotIn("mass accuracy", why)
            raws, fasta = _inputs(d)
            g, _ = _generate(d, gone, raws, fasta)
            self.assertNotEqual(g.returncode, 0)
            self.assertIn("cfg not found", g.stderr)
            self.assertNotIn("mass accuracy", g.stderr)
        with tempfile.TemporaryDirectory() as d:
            p = self._run(d, cfg=os.path.join(d, "nope.cfg"))
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("cfg not found", p.stderr)
            self.assertNotIn("mass accuracy", p.stderr)


    def test_with_no_existing_file_it_still_stops_the_chained_sbatch(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._run(d)
            self.assertEqual(p.returncode, run_search.SBATCH_NOT_WRITTEN, p.stderr)
            self.assertFalse(os.path.exists(os.path.join(d, "job.sh")))
            self.assertEqual(glob.glob(os.path.join(d, "job.sh.stale-*")), [])


class ProvenanceSaysWhatWasPassedTests(unittest.TestCase):
    """Round 2, item 4 (CLAUDE.md rule 1): every record describes exactly what DIA-NN is handed,
    and says "unverified" where its handling has not been measured."""

    def test_probe_path_mass_accuracy_is_reported_as_pinned(self):
        with tempfile.TemporaryDirectory() as d:
            raws, fasta = _inputs(d)
            p, _ = _generate(d, _cfg(d, PINNED_MA), raws, fasta)
            self.assertEqual(p.returncode, 0, p.stderr)
            info = json.loads(p.stdout)
            self.assertEqual(info["mass_acc"],
                             {"fixed": True, "ms1": 15.0, "ms2": 15.0,
                              "reason": "MS1 15.0 ppm / MS2 15.0 ppm, pinned in the cfg"})

    def test_override_with_a_pinned_window_says_pinned(self):
        """--allow-auto-mass-acc with `--window 7`: every step carries --window 7, and the record
        used to say "NOT pinned -- DIA-NN optimises per file"."""
        with tempfile.TemporaryDirectory() as d:
            raws, fasta = _inputs(d)
            p, out = _generate(d, _cfg(d, "--qvalue 0.01\n--window 7\n"), raws, fasta,
                               "--allow-auto-mass-acc")
            self.assertEqual(p.returncode, 0, p.stderr)
            info = json.loads(p.stdout)
            self.assertTrue(info["scan_window"]["source"].startswith("pinned in the cfg"),
                            info["scan_window"])
            self.assertEqual(info["scan_window"]["value"], 7)
            self.assertIn("--window 7", _read(os.path.join(out, "step2_firstpass.sbatch")))
            self.assertFalse(info["mass_acc"]["fixed"])

    def test_single_shot_records_an_unusable_window_as_passed_and_unverified(self):
        for body, passed in (("--window 7.0\n", ["--window 7.0"]),
                             ("--window 7\n--window 9\n", ["--window 7", "--window 9"])):
            with self.subTest(body=body), tempfile.TemporaryDirectory() as d:
                rec = run_search.scan_window_record("diann", _cfg(d, PINNED_MA + body), None)
                self.assertEqual(rec["passed"], passed)
                self.assertIn("unverified", rec["source"])
                self.assertNotIn("optimises it per run", rec["source"])
                self.assertIsNone(rec["value"])
        with tempfile.TemporaryDirectory() as d:
            rec = run_search.scan_window_record("diann", _cfg(d, PINNED_MA), None)
            self.assertIn("unverified", rec["source"])
            rec = run_search.scan_window_record("diann", _cfg(d, PINNED_MA + "--window 7\n", "p7.cfg"), None)
            self.assertEqual((rec["value"], rec["passed"]), (7, ["--window 7"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
