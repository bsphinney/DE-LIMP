#!/usr/bin/env python3
"""
Step 1b measures Orbitrap mass accuracy with DIA-NN, on the same representative runs and in the
same DIA-NN run per probe as the scan-window radius.

DIA-NN logs all three during calibration of a run whose --mass-acc, --mass-acc-ms1 and --window
are omitted. Verbatim from DIA-NN 2.7.0 on HIVE (srun job 23522741, 2026-09-16; Exploris 480
120k/15k run of the Set1-30-34-mouse cohort, full mouse predicted library, 32 threads):

    [0:36] Calibrating with mass accuracies 25 (MS1), 25 (MS2)
    [1:34] Scan window radius set to 7
    [1:35] Recommended MS1 mass accuracy setting: 4.1 ppm
    [2:27] Optimised mass accuracy: 14 ppm
    [2:45] Searching decoys
    ...
    Finished                                       (5:13 wall)

so a probe that stops at "Optimised mass accuracy" costs about half of the run's search.

The fake DIA-NN replays those captured logs (`<run>.log` beside each run) and behaves like the
real one where it matters, each point checked against DIA-NN 2.7.0 on HIVE:
  * with EITHER --mass-acc or --mass-acc-ms1 given it fixes both -- "WARNING: note the mass
    accuracy settings used by DIA-NN, automatic optimisation will not be performed as at least
    one of MS1/MS2 mass accuracies is user-provided", then "Mass accuracy will be fixed to ..."
    (srun 23528991: `--mass-acc-ms1 7` alone fixed MS2 at 2e-05) -- does not announce automatic
    optimisation and prints no "Optimised mass accuracy", but STILL prints "Recommended MS1 mass
    accuracy setting" (compare/pilot_override, pinned 20/7: "[1:50] Recommended MS1 mass
    accuracy setting: 4.3 ppm");
  * with --window N given it echoes "Scan window radius set to N" among the startup settings and
    infers none during calibration;
  * at "Searching decoys" it keeps running until the probe stops it.

Rebased onto step 1b's representative-run probe (fix/step1b-raw-dotnet-representative-probe): a
probe succeeds only when its run logged EVERYTHING it was asked to measure. A run that did not is
replaced by the next run nearest the median, exactly as a run with no radius is, and what it did
log is recorded but never pinned. A run whose flags FIX mass accuracy is not replaced: every other
run gets the same flags, so it stops at once, like DIA-NN's missing-.NET error.
"""
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import probe_window  # noqa: E402

GB = 1024 ** 3

# DIA-NN 2.7.0 log of run TT34, captured on HIVE (srun job 23522741) with mass accuracy and
# --window omitted: from the banner to the start of the search proper. Paths shortened; every
# other character is verbatim.
REAL_TT34 = """DIA-NN 2.7.0 Academia  (Data-Independent Acquisition by Neural Networks)
Compiled on Sep 15 2026 15:30:01
Current date and time: Wed Sep 16 19:07:09 2026
Logical CPU cores: 64
Thread number set to 32
Output will be filtered at 0.01 FDR
In silico digest will involve cuts at K*,R*
Min precursor m/z set to 357
Max precursor m/z set to 1105
Min precursor charge set to 2
Max precursor charge set to 4
Maximum number of missed cleavages set to 1
Min peptide length set to 7
Max peptide length set to 30
N-terminal methionine excision enabled
Cysteine carbamidomethylation enabled as a fixed modification
Maximum number of variable modifications set to 1
Modification UniMod:35 with mass delta 15.9949 at M will be considered as variable
DIA-NN will automatically optimise the mass accuracy for the first run of the experiment, use this mode for preliminary analyses only
WARNING: peptidoform scoring enabled because variable modifications have been declared; to disable, use --no-peptidoforms
The following variable modifications will be localised: UniMod:35

1 files will be processed
[0:00] Loading spectral library /scratch/lib/mouse.predicted.speclib
[0:03] Library annotated with sequence database(s): /data/fasta/UP000000589_10090/search.fasta
[0:04] Spectral library loaded: 22219 protein isoforms, 33078 protein groups and 3757675 precursors in 1783026 elution groups (targets and decoys).
[0:04] Loading protein annotations from FASTA /data/fasta/UP000000589_10090/search.fasta
[0:04] Annotating library proteins with information from the FASTA database
[0:04] Gene names missing for some isoforms
[0:04] Library contains 22193 proteins, and 21833 genes
[0:08] Initialising library

[0:17] File #1/1
[0:17] Loading run /data/raw/Ex01162023_12_TT34.raw
[0:35] Pre-processing...
[0:35] 1143 MS1 and 38842 MS2 scans in 1143 (inferred) and 1143 (encoded) cycles, 3757675 precursors in range
[0:36] Calibrating with mass accuracies 25 (MS1), 25 (MS2)
[1:34] RT window set to 1.11437
[1:34] Peak width: 3.436
[1:34] Scan window radius set to 7
[1:35] Recommended MS1 mass accuracy setting: 4.1 ppm
[2:27] Optimised mass accuracy: 14 ppm
[2:45] Searching decoys
"""

# The calibration lines of the other two representative runs of that cohort, same job, same
# settings, verbatim (the lines before them match TT34's apart from the run name and times).
# Across the three: radius 7, 7, 7; MS2 14, 17, 14 ppm; MS1 4.1, 4.3, 4.2 ppm.
REAL_CALIBRATION = {
    "Ex01162023_10_TT33": """[0:38] 1143 MS1 and 38858 MS2 scans in 1143 (inferred) and 1143 (encoded) cycles, 3757675 precursors in range
[0:39] Calibrating with mass accuracies 25 (MS1), 25 (MS2)
[1:49] RT window set to 1.13551
[1:49] Peak width: 3.444
[1:49] Scan window radius set to 7
[1:49] Recommended MS1 mass accuracy setting: 4.3 ppm
[2:48] Optimised mass accuracy: 17 ppm
[3:07] Searching decoys
""",
    "Ex01162023_8_TT32": """[0:37] 1144 MS1 and 38884 MS2 scans in 1144 (inferred) and 1144 (encoded) cycles, 3757675 precursors in range
[0:37] Calibrating with mass accuracies 25 (MS1), 25 (MS2)
[1:43] RT window set to 1.3045
[1:43] Peak width: 3.448
[1:43] Scan window radius set to 7
[1:44] Recommended MS1 mass accuracy setting: 4.2 ppm
[3:20] Optimised mass accuracy: 14 ppm
[3:50] Searching decoys
""",
}

# DIA-NN 2.7.0 with the flags pinned (--mass-acc 20 --mass-acc-ms1 7), same cohort, run TT33,
# from the step-1b branch's full-library timing run on HIVE (2026-09-16; library without Ox(M),
# hence 2,821,587 precursors): calibration still runs -- at its own 25 ppm -- but nothing is
# optimised, and the probe stopped it at the radius.
REAL_FIXED_LINES = """Mass accuracy will be fixed to 2e-05 (MS2) and 7e-06 (MS1)
[0:28] 1143 MS1 and 38858 MS2 scans in 1143 (inferred) and 1143 (encoded) cycles, 2821587 precursors in range
[0:28] Calibrating with mass accuracies 25 (MS1), 25 (MS2)
[1:39] RT window set to 1.14777
[1:39] Peak width: 3.388
[1:39] Scan window radius set to 7
"""


# DIA-NN 2.7.0 with ONLY --mass-acc-ms1 7 (and --window 7), run TT34, HIVE srun job 23528991
# (2026-09-16), verbatim from the settings to the search proper. One flag fixes both levels, MS2 at
# DIA-NN's 20 ppm default, and the MS1 recommendation is printed all the same.
REAL_ONE_FLAG_LINES = """Modification UniMod:35 with mass delta 15.9949 at M will be considered as variable
Scan window radius set to 7
WARNING: note the mass accuracy settings used by DIA-NN, automatic optimisation will not be performed as at least one of MS1/MS2 mass accuracies is user-provided
Mass accuracy will be fixed to 2e-05 (MS2) and 7e-06 (MS1)
WARNING: peptidoform scoring enabled because variable modifications have been declared; to disable, use --no-peptidoforms
The following variable modifications will be localised: UniMod:35

1 files will be processed
[0:17] Loading run /data/raw/Ex01162023_12_TT34.raw
[0:33] Calibrating with mass accuracies 25 (MS1), 25 (MS2)
[2:08] RT window set to 1.11437
[2:08] Recommended MS1 mass accuracy setting: 4.1 ppm
[2:35] Searching decoys
"""


def real_log(run):
    """The captured log of `run`, from the banner to the search proper."""
    if run == "Ex01162023_12_TT34":
        return REAL_TT34
    head = REAL_TT34[:REAL_TT34.index("[0:35] 1143 MS1")].replace("Ex01162023_12_TT34", run)
    return head + REAL_CALIBRATION[run]


FAKE_DIANN = r"""#!/bin/bash
[ -n "${FAKE_ARGV_LOG:-}" ] && echo "$@" >> "$FAKE_ARGV_LOG"
f=""; fixed=""; win=""
while [ $# -gt 0 ]; do
  case "$1" in
    --f) f="$2"; shift;;
    --mass-acc|--mass-acc-ms1) fixed=1; shift;;
    --window) win="$2"; shift;;
  esac; shift
done
case "$f" in
  *.raw|*.RAW)
    if [ -z "${DOTNET_ROOT:-}" ]; then
      echo "ERROR: cannot read .raw files, please download and install .NET Runtime 8: 8.0.17 or later https://dotnet.microsoft.com/en-us/download/dotnet/8.0 : 1"
      exit 0
    fi;;
esac
[ -n "$win" ] && echo "Scan window radius set to $win"
if [ -n "$fixed" ]; then
  if [ -n "${FAKE_FIXED_REWORDED:-}" ]; then
    # a DIA-NN build that words the fixed-accuracy notice differently
    echo "Mass tolerances pinned by the user"
  else
    echo "WARNING: note the mass accuracy settings used by DIA-NN, automatic optimisation will not be performed as at least one of MS1/MS2 mass accuracies is user-provided"
    echo "Mass accuracy will be fixed to 2e-05 (MS2) and 7e-06 (MS1)"
  fi
fi
while IFS= read -r line; do
  case "$line" in
    *"automatically optimise the mass accuracy"*|*"Optimised mass accuracy"*)
      [ -n "$fixed" ] && continue;;
    *"Scan window radius set to"*) [ -n "$win" ] && continue;;
    *"Searching decoys"*) echo "$line"; exec sleep "${FAKE_SEARCH_SLEEP:-20}";;
  esac
  case "$line" in *"Optimised mass accuracy"*) sleep "${FAKE_MS2_DELAY:-0}";; esac
  echo "$line"
done < "$f.log"
echo "Finished"
"""


def _exe(path, body):
    with open(path, "w") as fh:
        fh.write(body)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def _run(d, run, size, log=None):
    """A sparse .raw of `size` bytes with the DIA-NN log it replays beside it."""
    p = os.path.join(d, run + ".raw")
    with open(p, "wb") as fh:
        fh.truncate(size)
    with open(p + ".log", "w") as fh:
        fh.write(real_log(run) if log is None else log)
    return p


def _cohort(d):
    """The three representative runs by size: TT34 lower quartile, TT33 median, TT32 upper
    (their real sizes differ; only the order matters to the selection)."""
    return [_run(d, "Ex01162023_12_TT34", 13 * GB // 10), _run(d, "Ex01162023_10_TT33", 14 * GB // 10),
            _run(d, "Ex01162023_8_TT32", 15 * GB // 10)]


def no_ms2(log, radius=None):
    """`log` as a run that never reaches "Optimised mass accuracy" (and optionally logs another
    radius first). SYNTHETIC, derived from the captured logs."""
    out = re.sub(r"\[\d+:\d+\] Optimised mass accuracy: [0-9.]+ ppm\n", "", log)
    out = re.sub(r"\[\d+:\d+\] Searching decoys\n", "", out)
    if radius is not None:
        out = re.sub(r"Scan window radius set to \d+", f"Scan window radius set to {radius}", out)
    return out


def ms2_as(log, ppm):
    """`log` with its "Optimised mass accuracy" changed to `ppm`. SYNTHETIC."""
    return re.sub(r"Optimised mass accuracy: [0-9.]+ ppm", f"Optimised mass accuracy: {ppm} ppm", log)


def _fasta_lib(d):
    fasta = os.path.join(d, "db.fasta")
    with open(fasta, "w") as fh:
        fh.write(">sp|P1|X\nPEPTIDER\n")
    lib = os.path.join(d, "lib.speclib")
    with open(lib, "w") as fh:
        fh.write("lib")
    return fasta, lib


def _fake_dotnet_root(d):
    """A directory ensure_dotnet8.sh accepts as a .NET 8 >= 8.0.17 install (it also needs
    AspNetCore 8)."""
    root = os.path.join(d, "dotnet8")
    os.makedirs(root)
    _exe(os.path.join(root, "dotnet"),
         '#!/bin/bash\necho "Microsoft.AspNetCore.App 8.0.28 [%s/shared/Microsoft.AspNetCore.App]"\n'
         'echo "Microsoft.NETCore.App 8.0.28 [%s/shared/Microsoft.NETCore.App]"\n' % (root, root))
    return root


class LogParsingTests(unittest.TestCase):
    """read_log_line against DIA-NN 2.7.0's own words."""

    def _parse(self, text):
        found = {}
        for line in text.splitlines():
            for k, v in probe_window.read_log_line(line).items():
                found.setdefault(k, v)
        return found

    def test_the_real_auto_mode_log_yields_the_radius_and_both_mass_accuracies(self):
        self.assertEqual(self._parse(REAL_TT34), {"radius": 7, "ms1_ppm": 4.1, "ms2_ppm": 14.0})

    def test_a_pinned_run_reports_no_mass_accuracy(self):
        """'Calibrating with mass accuracies 25 (MS1), 25 (MS2)' is DIA-NN's calibration-stage
        tolerance, printed whether or not anything is optimised -- never a measurement."""
        found = self._parse(REAL_FIXED_LINES)
        self.assertEqual(found, {"radius": 7})
        self.assertTrue(probe_window.FIXED_ACC_RE.search(REAL_FIXED_LINES))

    def _read(self, text, measure):
        r = probe_window.LogReader(measure)
        for line in text.splitlines():
            r.feed(line)
        return r

    def test_a_pinned_runs_ms1_recommendation_is_not_a_measurement(self):
        """DIA-NN 2.7.0 prints "Recommended MS1 mass accuracy setting" with the flags pinned too.
        The reader accepts mass accuracy only from a run that announced automatic optimisation,
        so a pinned run yields none -- and is flagged at once, not after the search."""
        r = self._read(REAL_ONE_FLAG_LINES, ["mass-acc"])
        self.assertNotIn("ms1_ppm", r.found)
        self.assertTrue(r.found.get("fixed_mass_acc"))
        self.assertTrue(r.done())
        auto = self._read(REAL_TT34, ["window", "mass-acc"])
        self.assertEqual({k: auto.found[k] for k in ("radius", "ms1_ppm", "ms2_ppm")},
                         {"radius": 7, "ms1_ppm": 4.1, "ms2_ppm": 14.0})
        self.assertFalse(auto.found.get("fixed_mass_acc"))

    def test_a_reworded_fixed_notice_still_cannot_pass_as_auto_mode(self):
        """Positive confirmation, not a blacklist: if a DIA-NN build rewords "Mass accuracy will be
        fixed to", the end of the settings block without "automatically optimise the mass
        accuracy" still says the run is not measuring anything.

        The two are recorded SEPARATELY. "DIA-NN said it is fixing it" is fixed by removing the
        flags; "the settings block ended without the announcement" may instead mean the wording
        changed between DIA-NN releases, and its message must not send the reader after flags
        that are already correct."""
        reworded = REAL_ONE_FLAG_LINES.replace(
            "Mass accuracy will be fixed to 2e-05 (MS2) and 7e-06 (MS1)", "Tolerances pinned") \
            .replace("WARNING: note the mass accuracy settings used by DIA-NN, automatic "
                     "optimisation will not be performed as at least one of MS1/MS2 mass "
                     "accuracies is user-provided\n", "")
        self.assertFalse(probe_window.FIXED_ACC_RE.search(reworded))
        r = self._read(reworded, ["mass-acc"])
        self.assertNotIn("ms1_ppm", r.found)
        self.assertTrue(r.not_optimising())
        self.assertTrue(r.done())
        self.assertTrue(r.found.get("no_auto_announcement"))
        self.assertFalse(r.found.get("fixed_mass_acc"),
                         "a reworded notice was reported as DIA-NN saying it fixed the tolerance")
        # ...and the notice DIA-NN 2.7.0 really prints is the other one, never this one
        said = self._read(REAL_ONE_FLAG_LINES, ["mass-acc"])
        self.assertTrue(said.found.get("fixed_mass_acc"))
        self.assertFalse(said.found.get("no_auto_announcement"))

    def test_the_measurement_is_the_median_of_each_level_as_printed(self):
        """No rounding: the MEASURED number is one DIA-NN printed for a run. What is pinned may
        then be the SOP floor instead -- that is the next test; this one is about the median."""
        probes = [{"ms2_ppm": 14.0, "ms1_ppm": 4.1}, {"ms2_ppm": 13.0, "ms1_ppm": 4.3},
                  {"ms2_ppm": 15.0, "ms1_ppm": 3.9}]
        pin = probe_window.pin_mass_acc(probes)
        self.assertEqual((pin["measured_ms2_ppm"], pin["measured_ms1_ppm"]), (14.0, 4.1))
        self.assertEqual(pin["ms2_per_run"], [14.0, 13.0, 15.0])
        self.assertFalse(pin["agree"])
        one = probe_window.pin_mass_acc([{"ms2_ppm": 14.0, "ms1_ppm": 4.1}])
        self.assertEqual((one["measured_ms2_ppm"], one["measured_ms1_ppm"]), (14.0, 4.1))

    def test_the_median_absorbs_a_run_whose_optimisation_moved(self):
        """Measured on HIVE, same three runs, same library, radius 7 both times: with --window
        auto, TT32 gave MS2 14 / MS1 4.2; with --window 7 pinned it gave MS2 12 / MS1 4.5 (TT34
        and TT33 were unchanged). That is not noise: the pinned-window values reproduced exactly
        in two separate sruns (23524590 and 23526415). Pinning the window changes the conditions
        DIA-NN optimises under. The pinned MS2 does not move; a measured MS1 moves by 0.1 ppm."""
        auto = probe_window.pin_mass_acc([{"ms2_ppm": 14.0, "ms1_ppm": 4.1},
                                          {"ms2_ppm": 17.0, "ms1_ppm": 4.3},
                                          {"ms2_ppm": 14.0, "ms1_ppm": 4.2}])
        pinned_window = probe_window.pin_mass_acc([{"ms2_ppm": 14.0, "ms1_ppm": 4.1},
                                                   {"ms2_ppm": 17.0, "ms1_ppm": 4.3},
                                                   {"ms2_ppm": 12.0, "ms1_ppm": 4.5}])
        self.assertEqual((auto["measured_ms2_ppm"], auto["measured_ms1_ppm"]), (14.0, 4.2))
        self.assertEqual((pinned_window["measured_ms2_ppm"],
                          pinned_window["measured_ms1_ppm"]), (14.0, 4.3))
        # both are tighter than the SOP, so both SEARCH at the SOP -- the measurement is what
        # moved, and it is the measurement this test is about
        self.assertEqual(auto["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")
        self.assertEqual(pinned_window["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")

    def test_a_documented_level_is_pinned_as_documented_and_still_recorded(self):
        """120k MS1 has a README tier (7 ppm). The measured MS1 (4.1-4.5) is kept as evidence but
        not pinned: pinned, DIA-NN 2.7.0 warned on every pass that it "deviates significantly
        from the value recommended (7 ppm) for the Orbitrap resolution of this run (120000)"."""
        probes = [{"ms2_ppm": 14.0, "ms1_ppm": 4.1}, {"ms2_ppm": 17.0, "ms1_ppm": 4.3},
                  {"ms2_ppm": 12.0, "ms1_ppm": 4.5}]
        pin = probe_window.pin_mass_acc(probes, documented={"ms1_ppm": 7.0})
        self.assertEqual(pin["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")
        self.assertEqual(pin["ms1_per_run"], [4.1, 4.3, 4.5])
        # the documented level is pinned as given and NEVER floored; the measurement it did not
        # use is still recorded
        self.assertFalse(pin["floored"]["ms1_ppm"])
        self.assertEqual(pin["measured_ms1_ppm"], 4.3)
        self.assertEqual(pin["pinned_ms1_ppm"], 7.0)
        self.assertIn("documented", pin["sources"]["ms1_ppm"])
        self.assertIn("4.3 ppm (recorded, not used)", pin["sources"]["ms1_ppm"])
        self.assertIn("median", pin["sources"]["ms2_ppm"])
        self.assertFalse(pin["agree"], "agreement is over the MEASURED levels")
        # a documented level needs no per-run value at all
        pin = probe_window.pin_mass_acc([{"ms2_ppm": 14.0}], documented={"ms1_ppm": 7.0})
        self.assertEqual(pin["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")
        self.assertIsNone(pin["measured_ms1_ppm"])
        self.assertIsNone(probe_window.pin_mass_acc([{"ms1_ppm": 4.1}], documented={"ms1_ppm": 7}))

    def test_no_pin_unless_every_run_logged_both_levels(self):
        self.assertIsNone(probe_window.pin_mass_acc(
            [{"ms2_ppm": 14.0, "ms1_ppm": 4.1}, {"ms2_ppm": None, "ms1_ppm": 4.3}]))
        self.assertIsNone(probe_window.pin_mass_acc([]))


class ProbeMassAccTests(unittest.TestCase):
    """probe_window.py end to end, replaying the captured logs."""

    def _probe(self, d, raws, measure=("window", "mass-acc"), extra="", env_extra=None,
               timeout=60, more=(), after=()):
        fasta, lib = _fasta_lib(d)
        diann = _exe(os.path.join(d, "diann"), FAKE_DIANN)
        cfg = os.path.join(d, "resolved.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01")                    # no trailing newline, on purpose
        argv_log = os.path.join(d, "argv.txt")
        env = dict(os.environ, DOTNET_ROOT="/opt/fake-dotnet", FAKE_ARGV_LOG=argv_log)
        env.update(env_extra or {})
        argv = [sys.executable, os.path.join(SCRIPTS, "probe_window.py"), "--diann", diann,
                "--raw", *raws, "--fasta", fasta, "--lib", lib, "--threads", "8",
                "--timeout", str(timeout), "--write-cfg", cfg,
                "--workdir", os.path.join(d, "w"), "--measure", *measure, *more]
        if extra:
            argv += ["--extra", extra]
        if after:                                        # DIA-NN flags as the chain passes them
            argv += ["--", *after]
        t0 = time.time()
        p = subprocess.run(argv, capture_output=True, text=True, env=env, timeout=240)
        return p, cfg, argv_log, time.time() - t0

    def test_one_diann_per_run_measures_window_and_mass_accuracy_and_stops_there(self):
        with tempfile.TemporaryDirectory() as d:
            p, cfg, argv_log, secs = self._probe(d, _cohort(d))
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(out["measured"], ["window", "mass-acc"])
            self.assertEqual(out["stopped_because"], "measured")
            self.assertFalse(out["incomplete"])
            # every run's values are in the job log, not only in the JSON
            self.assertRegex(p.stderr, r"Ex01162023_10_TT33\.raw.*radius 7.*MS2 17 ppm.*MS1 4\.3 ppm")
            per_run = {os.path.basename(x["file"]): (x["radius"], x["ms2_ppm"], x["ms1_ppm"])
                       for x in out["probes"]}
            self.assertEqual(per_run["Ex01162023_12_TT34.raw"], (7, 14.0, 4.1))
            ma = out["mass_acc"]
            self.assertEqual(ma["pin_as"], "--mass-acc %s --mass-acc-ms1 %s" % (
                probe_window.ppm_text(ma["ms2_ppm"]), probe_window.ppm_text(ma["ms1_ppm"])))
            self.assertEqual(len(open(argv_log).read().splitlines()), 3,
                             "one DIA-NN per run, radius and mass accuracy from the same one")
            # each fake keeps "searching" for 20 s once it reaches Searching decoys
            self.assertLess(secs, 30, "a probe waited for the search instead of stopping")
            txt = open(cfg).read()
            self.assertIn("--qvalue 0.01\n", txt)
            for flag in ("--window", "--mass-acc", "--mass-acc-ms1"):
                self.assertEqual(len(re.findall(r"^%s " % re.escape(flag), txt, re.M)), 1, txt)
            self.assertIn(ma["pin_as"].split(" --mass-acc-ms1 ")[0] + "\n", txt)

    def test_the_probe_keeps_reading_after_the_radius_until_the_ms2_line(self):
        """On the real log the MS2 value comes ~50 s after the radius. A window-only probe must
        not wait for it; a mass-accuracy probe must."""
        with tempfile.TemporaryDirectory() as d:
            raws = [_run(d, "Ex01162023_12_TT34", GB)]
            p, _, _, secs = self._probe(d, raws, measure=("window",),
                                        env_extra={"FAKE_MS2_DELAY": "6"})
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertLess(secs, 5, "a window-only probe waited for mass accuracy")
            self.assertIsNone(json.loads(p.stdout)["mass_acc"])
            # mass accuracy needs MASS_ACC_MIN_RUNS runs, so this half gets the cohort
            p, _, _, secs = self._probe(d, _cohort(d), env_extra={"FAKE_MS2_DELAY": "3"})
            self.assertEqual(p.returncode, 0, p.stderr)
            ma = json.loads(p.stdout)["mass_acc"]
            self.assertEqual(ma["measured_ms2_ppm"], 14.0)
            self.assertEqual(ma["pinned_ms2_ppm"], 20.0)

    def test_a_run_that_logs_no_mass_accuracy_is_replaced_and_what_it_did_log_is_not_pinned(self):
        """Step 1b's rule for a run with no radius, applied to a run with no mass accuracy: the
        next run nearest the median takes its place. The failed run DID log a radius (9 here);
        a run counts only when it logged everything asked, so that 9 is recorded, not pinned --
        otherwise the radius and the mass accuracy would describe different sets of runs."""
        with tempfile.TemporaryDirectory() as d:
            raws = _cohort(d) + [_run(d, "Ex01162023_14_TT31", 16 * GB // 10,
                                      log=ms2_as(real_log("Ex01162023_12_TT34"), 15)
                                      .replace("Ex01162023_12_TT34", "Ex01162023_14_TT31"))]
            # by size: TT34, TT33, TT32, TT31 -> median TT32, quartiles TT34 and TT31, TT33 reserve
            with open(raws[2] + ".log", "w") as fh:
                fh.write(no_ms2(real_log("Ex01162023_8_TT32"), radius=9))
            p, cfg, _, _ = self._probe(d, raws)
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            names = [os.path.basename(x["file"]) for x in out["probes"]]
            self.assertEqual(names, ["Ex01162023_8_TT32.raw", "Ex01162023_12_TT34.raw",
                                     "Ex01162023_14_TT31.raw", "Ex01162023_10_TT33.raw"])
            self.assertEqual(out["failed"], ["Ex01162023_8_TT32.raw"])
            self.assertEqual(out["probes"][0]["radius"], 9, "what the failed run logged is recorded")
            self.assertEqual(out["probes"][0]["missing"], ["ms2_ppm"])
            self.assertEqual(out["probes"][3]["role"], "reserve")
            self.assertEqual(out["radii"], [7, 7, 7], "the failed run's radius was used")
            self.assertEqual(out["mass_acc"]["ms2_per_run"], [14.0, 15.0, 17.0])
            self.assertEqual(out["mass_acc"]["measured_ms2_ppm"], 15.0)
            self.assertEqual(out["mass_acc"]["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")
            self.assertFalse(out["incomplete"])
            self.assertIn("--window 7\n", open(cfg).read())

    def test_when_no_run_logs_mass_accuracy_it_fails_names_them_and_pins_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            raws = _cohort(d)
            for r in raws:
                with open(r + ".log", "w") as fh:
                    fh.write(no_ms2(real_log(os.path.basename(r)[:-4])))
            p, cfg, _, _ = self._probe(d, raws)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("Could not read the scan-window radius and mass accuracy", p.stderr)
            self.assertIn("Ex01162023_10_TT33.raw", p.stderr)
            self.assertIn("'Optimised mass accuracy' (MS2)", p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(out["stopped_because"], "max_failures")
            self.assertIsNone(out["mass_acc"])
            self.assertIsNone(out["window_radius"], "a radius was pinned from a failed step")
            self.assertEqual([x["radius"] for x in out["probes"]], [7, 7, 7])
            self.assertNotIn("--mass-acc", open(cfg).read())
            self.assertNotIn("--window", open(cfg).read())

    def test_a_documented_ms1_is_pinned_as_given(self):
        with tempfile.TemporaryDirectory() as d:
            p, cfg, _, _ = self._probe(d, _cohort(d), measure=("mass-acc",),
                                       more=["--ms1-ppm", "7"])
            self.assertEqual(p.returncode, 0, p.stderr)
            ma = json.loads(p.stdout)["mass_acc"]
            # MS2 per run 14 / 17 / 14 -> measured 14, tighter than the SOP so pinned at 20;
            # MS1 documented, pinned as given and never floored
            self.assertEqual(ma["measured_ms2_ppm"], 14.0)
            self.assertEqual(ma["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")
            self.assertEqual(ma["floored"], {"ms2_ppm": True, "ms1_ppm": False})
            self.assertIn("measured 14 ppm", ma["sources"]["ms2_ppm"])
            self.assertIn("pinned at 20 ppm", ma["sources"]["ms2_ppm"])
            self.assertRegex(p.stderr, r"MS2: measured 14 ppm, PINNED 20 ppm")
            # probing order: the median run (TT33) first, then TT34 and TT32
            self.assertEqual(ma["ms1_per_run"], [4.3, 4.1, 4.2])
            self.assertIn("--mass-acc-ms1 7\n", open(cfg).read())

    def test_documenting_both_levels_leaves_nothing_to_measure(self):
        with tempfile.TemporaryDirectory() as d:
            p, _, argv_log, _ = self._probe(d, _cohort(d), measure=("mass-acc",),
                                            more=["--ms1-ppm", "7", "--ms2-ppm", "15"])
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("nothing to measure", p.stderr)
            self.assertFalse(os.path.exists(argv_log), "DIA-NN ran with nothing to measure")

    def test_a_run_that_never_announced_auto_mode_yields_no_mass_accuracy(self):
        """The regression the fake used to hide: real DIA-NN prints the MS1 recommendation with
        the flags pinned. If the fixed-accuracy notice were ever reworded, the probe must still
        not record that recommendation as a measurement, and must stop at once -- and it must say
        which of the two faults it saw, because they have different fixes."""
        with tempfile.TemporaryDirectory() as d:
            raws = [_run(d, "Ex01162023_12_TT34", GB)]
            p, cfg, _, secs = self._probe(d, raws, measure=("mass-acc",),
                                          extra="--mass-acc 20 --mass-acc-ms1 7", timeout=120,
                                          env_extra={"FAKE_FIXED_REWORDED": "1"})
            self.assertNotEqual(p.returncode, 0)
            rec = json.loads(p.stdout)["probes"][0]
            self.assertIsNone(rec["ms1_ppm"], "a pinned run's recommendation became the pin")
            self.assertTrue(rec["mass_acc_no_auto_announcement"])
            self.assertFalse(rec["mass_acc_fixed_by_flags"],
                             "a reworded notice was reported as DIA-NN saying it fixed it")
            # the flags here DO carry both, so the message names them rather than the wording
            self.assertIn("--mass-acc and --mass-acc-ms1", p.stderr)
            self.assertLess(secs, 15, "the probe waited for a search that measures nothing")
            self.assertNotIn("--mass-acc", open(cfg).read())

    def test_a_reordered_settings_block_does_not_blame_flags_that_are_correct(self):
        """A DIA-NN release that rewords or reorders the auto-optimisation announcement lands in
        the same place as a pinned run: the settings block ends without it. But nothing is wrong
        with the flags -- there is no --mass-acc or --mass-acc-ms1 to remove -- so the message
        must say so and point at the log wording instead. Otherwise the one instruction the
        reader gets is to remove flags that are not there."""
        with tempfile.TemporaryDirectory() as d:
            raws = _cohort(d)
            for r in raws:                                  # the announcement, reworded
                text = open(r + ".log").read().replace(
                    "DIA-NN will automatically optimise the mass accuracy for the first run of "
                    "the experiment", "DIA-NN will auto-tune the mass tolerance for run 1")
                with open(r + ".log", "w") as fh:
                    fh.write(text)
            p, cfg, argv_log, secs = self._probe(d, raws, measure=("mass-acc",), timeout=120)
            self.assertNotEqual(p.returncode, 0)
            out = json.loads(p.stdout)
            self.assertEqual(out["stopped_because"], "environment")
            self.assertEqual(len(out["probes"]), 1, "a run was replaced over the log wording")
            self.assertTrue(out["probes"][0]["mass_acc_no_auto_announcement"])
            self.assertFalse(out["probes"][0]["mass_acc_fixed_by_flags"])
            self.assertIn("NEITHER --mass-acc NOR --mass-acc-ms1", p.stderr)
            self.assertIn("AUTO_ACC_RE", p.stderr)
            self.assertNotIn("remove them", p.stderr)
            self.assertIsNone(out["mass_acc"])
            self.assertNotIn("--mass-acc", open(cfg).read())

    def test_flags_that_fix_mass_accuracy_fail_fast_with_the_reason(self):
        """DIA-NN optimises only what is omitted; a probe must not sit out its whole budget
        waiting for lines a pinned run never prints. And it is not a property of the run --
        every other run gets the same flags -- so nothing is replaced: it stops at once, as for
        DIA-NN's missing-.NET error. The flags may come as --extra or after `--` (the chain)."""
        for how in ("extra", "after"):
            with self.subTest(how), tempfile.TemporaryDirectory() as d:
                flags = ["--mass-acc", "20", "--mass-acc-ms1", "7"]
                p, _, argv_log, secs = self._probe(
                    d, _cohort(d), timeout=120, **({"extra": " ".join(flags)} if how == "extra"
                                                   else {"after": flags}))
                self.assertNotEqual(p.returncode, 0)
                self.assertLess(secs, 15)
                self.assertIn("fix", p.stderr)
                out = json.loads(p.stdout)
                self.assertEqual(len(out["probes"]), 1, "a run was replaced over the flags")
                self.assertEqual(out["stopped_because"], "environment")
                self.assertTrue(out["probes"][0]["mass_acc_fixed_by_flags"])
                self.assertIn("--mass-acc 20 --mass-acc-ms1 7", open(argv_log).read())


class PinnedWindowGuardTests(unittest.TestCase):
    """With --window N given, DIA-NN 2.7.0 echoes "Scan window radius set to N" at startup (HIVE,
    2026-09-16). A probe asked to measure the window with --window in --extra would read that
    echo back as a measurement -- the cfg's number, reported as DIA-NN's."""

    def test_measuring_a_window_that_the_flags_pin_is_refused(self):
        """--extra, or after `--` as step 1b passes the cfg's flags."""
        for tail in (["--extra", "--window 9"], ["--", "--qvalue", "0.01", "--window", "9"]):
            with self.subTest(tail), tempfile.TemporaryDirectory() as d:
                fasta, lib = _fasta_lib(d)
                raw = _run(d, "Ex01162023_12_TT34", GB)
                argv_log = os.path.join(d, "argv.txt")
                p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "probe_window.py"),
                                    "--diann", _exe(os.path.join(d, "diann"), FAKE_DIANN),
                                    "--raw", raw, "--fasta", fasta, "--lib", lib,
                                    "--workdir", os.path.join(d, "w"), *tail],
                                   capture_output=True, text=True, timeout=60,
                                   env=dict(os.environ, DOTNET_ROOT="/opt/fake-dotnet",
                                            FAKE_ARGV_LOG=argv_log))
                self.assertNotEqual(p.returncode, 0)
                self.assertIn("--window", p.stderr)
                self.assertNotIn('"window_radius": 9', p.stdout)
                self.assertFalse(os.path.exists(argv_log), "DIA-NN ran to echo a pinned window")

    def test_measuring_only_mass_accuracy_under_a_pinned_window_is_allowed(self):
        with tempfile.TemporaryDirectory() as d:
            fasta, lib = _fasta_lib(d)
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "probe_window.py"),
                                "--diann", _exe(os.path.join(d, "diann"), FAKE_DIANN),
                                "--raw", *_cohort(d), "--fasta", fasta, "--lib", lib,
                                "--workdir", os.path.join(d, "w"), "--measure", "mass-acc",
                                "--ms1-ppm", "7", "--", "--window", "7"],
                               capture_output=True, text=True, timeout=120,
                               env=dict(os.environ, DOTNET_ROOT="/opt/fake-dotnet"))
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertIsNone(out["window_radius"])
            self.assertEqual(out["radii"], [])
            self.assertEqual(out["mass_acc"]["measured_ms2_ppm"], 14.0)
            self.assertEqual(out["mass_acc"]["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")


class Step1bMassAccChainTests(unittest.TestCase):
    """diann_parallel.py generates step 1b and steps 2-5 for a cfg estimate_params.py planned as
    `measure_with_diann`; step 1b is then executed with bash, as a compute node would."""

    def _chain(self, d, raws, cfg_extra=""):
        fasta, _ = _fasta_lib(d)
        cfg = os.path.join(d, "params.cfg")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "diann", "--acquisition", "DIA",
                            "--instrument", "Orbitrap Exploris 480",
                            "--ms1-resolution", "120000", "--ms2-resolution", "15000",
                            "--precursor-mz-range", "357", "1105", "--out", cfg],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        if cfg_extra:
            with open(cfg, "a") as fh:
                fh.write(cfg_extra)
        diann = _exe(os.path.join(d, "diann"), FAKE_DIANN)
        out = os.path.join(d, "out")
        env = {k: v for k, v in os.environ.items() if k != "DOTNET_ROOT"}
        env["PROTEOMICS_DOTNET_DIR"] = _fake_dotnet_root(d)
        g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                            "--diann", diann, "--raw", *raws, "--fasta", fasta,
                            "--out", out, "--cfg", cfg, "--threads-per-file", "8"],
                           capture_output=True, text=True, env=env, timeout=120)
        self.assertEqual(g.returncode, 0, g.stderr)
        open(os.path.join(out, "step1.predicted.speclib"), "w").write("lib")
        return out, json.loads(g.stdout)

    def _step(self, out, name):
        return open(os.path.join(out, name)).read()

    def _run_step1b(self, d, out):
        env = {k: v for k, v in os.environ.items()
               if k not in ("DOTNET_ROOT", "PROTEOMICS_DOTNET_DIR")}
        env["FAKE_ARGV_LOG"] = os.path.join(d, "argv.txt")
        return subprocess.run(["bash", os.path.join(out, "step1b_window.sbatch")], cwd=out,
                              capture_output=True, text=True, env=env, timeout=240)

    def test_step1b_measures_mass_accuracy_and_every_step_reads_it(self):
        with tempfile.TemporaryDirectory() as d:
            out, info = self._chain(d, _cohort(d))
            self.assertEqual(info["step1b_measures"], ["window", "mass-acc"])
            self.assertEqual(info["parallel_safe"]["code"], "probe")
            massacc = os.path.join(out, "massacc.txt")
            # provenance says what will be passed: measured at run time, the documented MS1 known
            # now -- not upstream's record for an omitted flag ("DIA-NN calibrates it itself")
            ma = info["mass_acc"]
            self.assertEqual(ma["value_file"], massacc)
            self.assertEqual(ma["evidence_file"], os.path.join(out, "window.json"))
            self.assertTrue(ma["measured"])
            # "measured: true" means we MEASURE it, not that the search runs at the measured
            # number: a measured level is floored at the SOP, and the record has to say so
            self.assertEqual(ma["sop_floor"], {"--mass-acc": 20.0, "--mass-acc-ms1": 7.0})
            self.assertIn("max(measured, SOP)", ma["floor_note"])
            self.assertIn("measured_ms2_ppm", ma["floor_note"])
            self.assertIn("pinned_ms2_ppm", ma["floor_note"])
            self.assertEqual((ma["ms1"], ma["ms2"]), (7, None))
            self.assertEqual(ma["documented"], {"--mass-acc-ms1": 7})
            self.assertIn("step 1b", ma["source"])
            self.assertNotIn("calibrates it itself", json.dumps(ma))
            self.assertEqual(info["resolved_params"]["produced"], "runtime")
            s1b = self._step(out, "step1b_window.sbatch")
            self.assertIn("--measure window mass-acc", s1b)
            self.assertIn("--ms1-ppm 7", s1b, "120k MS1 has a README tier; step 1b must keep it")
            self.assertNotIn("--ms2-ppm", s1b)
            for name in ("step2_firstpass.sbatch", "step3_assembly.sbatch",
                         "step4_finalpass.sbatch", "step5_report.sbatch"):
                body = self._step(out, name)
                self.assertIn(f"$(cat {massacc})", body, name)
                self.assertIn(f"--window $(cat {os.path.join(out, 'window.txt')})", body, name)
                self.assertNotRegex(body, r"--mass-acc(-ms1)? \d", name)
            self.assertNotIn("--mass-acc", self._step(out, "step1_libpred.sbatch"))

            p = self._run_step1b(d, out)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            w = json.load(open(os.path.join(out, "window.json")))
            pin = w["mass_acc"]["pin_as"]
            self.assertEqual(w["mass_acc"]["measured_ms2_ppm"], 14.0,
                             "the measurement must survive the SOP floor in the evidence")
            # MS2: the median of the three real runs (14/17/14); MS1: the README's 120k tier,
            # not the measured 4.1/4.3/4.2 -- which DIA-NN 2.7.0 warns about when pinned
            self.assertEqual(pin, "--mass-acc 20 --mass-acc-ms1 7")
            # probing order: the median run (TT33) first, then TT34 and TT32
            self.assertEqual(w["mass_acc"]["ms1_per_run"], [4.3, 4.1, 4.2])
            self.assertEqual(open(os.path.join(out, "window.txt")).read().strip(), "7")
            self.assertEqual(open(massacc).read().strip(), pin)
            self.assertEqual(len(pin.split()), 4, "massacc.txt must expand to two flags")
            self.assertIn("mass accuracy = " + pin, p.stdout)
            resolved = open(os.path.join(out, "params.resolved.cfg")).read()
            for flag in ("--mass-acc", "--mass-acc-ms1", "--window"):
                self.assertEqual(len(re.findall(r"^%s " % re.escape(flag), resolved, re.M)), 1,
                                 resolved)
            # the probes ran DIA-NN in auto mode: nothing pinned mass accuracy on them
            for call in open(os.path.join(d, "argv.txt")).read().splitlines():
                self.assertNotIn("--mass-acc", call)

    def test_a_pinned_window_leaves_step1b_measuring_mass_accuracy_only(self):
        with tempfile.TemporaryDirectory() as d:
            out, info = self._chain(d, _cohort(d), cfg_extra="--window 7\n")
            self.assertEqual(info["step1b_measures"], ["mass-acc"])
            # the window is the cfg's, and provenance says so -- not "measured by step 1b"
            self.assertTrue(info["scan_window"]["source"].startswith("pinned in the cfg"),
                            info["scan_window"])
            self.assertIn("--measure mass-acc ", self._step(out, "step1b_window.sbatch"))
            s2 = self._step(out, "step2_firstpass.sbatch")
            self.assertIn("--window 7", s2)
            self.assertNotIn("window.txt", s2)
            self.assertIn("massacc.txt", s2)
            p = self._run_step1b(d, out)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertFalse(os.path.exists(os.path.join(out, "window.txt")))
            # the probes ran under the cfg's --window, as steps 2-5 will
            for call in open(os.path.join(d, "argv.txt")).read().splitlines():
                self.assertIn("--window 7", call)
            # DIA-NN echoes the pinned --window at startup; that is the cfg, not a measurement
            w = json.load(open(os.path.join(out, "window.json")))
            self.assertEqual([x["radius"] for x in w["probes"]], [None, None, None])
            self.assertIsNone(w["window_radius"])
            self.assertTrue(os.path.getsize(os.path.join(out, "massacc.txt")) > 0)
            resolved = open(os.path.join(out, "params.resolved.cfg")).read()
            self.assertEqual(resolved.count("--window 7"), 1, resolved)

    def _run_step(self, d, out, name, task=0):
        env = {k: v for k, v in os.environ.items()
               if k not in ("DOTNET_ROOT", "PROTEOMICS_DOTNET_DIR")}
        argv_log = os.path.join(d, "argv_%s.txt" % name.split("_")[0])
        env.update(FAKE_ARGV_LOG=argv_log, SLURM_ARRAY_TASK_ID=str(task), FAKE_SEARCH_SLEEP="0")
        p = subprocess.run(["bash", os.path.join(out, name)], cwd=out, capture_output=True,
                           text=True, env=env, timeout=120)
        return p, argv_log

    STEPS = ("step2_firstpass.sbatch", "step3_assembly.sbatch", "step4_finalpass.sbatch",
             "step5_report.sbatch")

    def test_steps_2_to_5_refuse_to_run_without_what_step1b_measured(self):
        """`$(cat massacc.txt)` expands to NOTHING when the file is missing, and DIA-NN then runs
        in auto mode, optimising per file -- what the chain exists to prevent -- and still writes
        a .quant, so must_exist passes. Reproduced in review: step 1b failed (after its own
        `rm -f massacc.txt`), references/watcher.md says to resubmit downstream steps, and
        `sbatch step2_firstpass.sbatch` ran DIA-NN with no mass-accuracy flag at all."""
        for cfg_extra, missing in (("", ("window.txt", "massacc.txt")),     # checked in this order
                                   ("--window 7\n", ("massacc.txt",))):
            with tempfile.TemporaryDirectory() as d:
                out, _ = self._chain(d, _cohort(d), cfg_extra=cfg_extra)
                for name in self.STEPS:
                    p, argv_log = self._run_step(d, out, name)
                    self.assertNotEqual(p.returncode, 0, name)
                    self.assertIn("step 1b", p.stderr, name)
                    self.assertIn(missing[0], p.stderr, name)
                    self.assertFalse(os.path.exists(argv_log),
                                     f"{name} ran DIA-NN without the values step 1b measures")
                # present but not a measurement -- e.g. step 1b printed "None" into it
                with open(os.path.join(out, "massacc.txt"), "w") as fh:
                    fh.write("None\n")
                with open(os.path.join(out, "window.txt"), "w") as fh:
                    fh.write("7\n")
                p, argv_log = self._run_step(d, out, "step2_firstpass.sbatch")
                self.assertNotEqual(p.returncode, 0)
                self.assertIn("massacc.txt does not hold", p.stderr)
                self.assertFalse(os.path.exists(argv_log), "massacc.txt 'None' reached DIA-NN")

    def test_steps_2_to_5_refuse_a_massacc_file_that_is_not_a_plausible_tolerance(self):
        """MEASURED_FILE_RE checks the SHAPE of massacc.txt -- two flags, two numbers -- and a
        shape check cannot tell 14 ppm from 999999: `--mass-acc 999999 --mass-acc-ms1 0.001`
        passed every guard in the chain. The probe can no longer write such a line, so this is
        the guard for every other way the file can hold one (hand-edited, copied from another
        cohort, left by an older version). Steps 2-5 splice `$(cat massacc.txt)` straight onto a
        DIA-NN command line and nothing downstream looks at the number again."""
        for line in ("--mass-acc 999999 --mass-acc-ms1 0.001",
                     "--mass-acc 0.4 --mass-acc-ms1 7",
                     "--mass-acc 60 --mass-acc-ms1 7",
                     "--mass-acc 14 --mass-acc-ms1 40"):
            with self.subTest(line), tempfile.TemporaryDirectory() as d:
                out, _ = self._chain(d, _cohort(d))
                with open(os.path.join(out, "window.txt"), "w") as fh:
                    fh.write("7\n")
                with open(os.path.join(out, "massacc.txt"), "w") as fh:
                    fh.write(line + "\n")
                for name in self.STEPS:
                    p, argv_log = self._run_step(d, out, name)
                    self.assertNotEqual(p.returncode, 0, name)
                    self.assertIn("outside the plausible band", p.stderr, name)
                    self.assertFalse(os.path.exists(argv_log),
                                     f"{name} searched at {line}")

    def test_a_measured_pin_still_passes_the_band_guard(self):
        """The guard above must not refuse what step 1b really measures, nor the values a maintainer
        may legitimately pin instead (the open question: measured 14 vs DIA-NN's own 25 for a 15k
        MS2, vs the pilot's hand-set 20). All of them are inside the band."""
        for line in ("--mass-acc 14 --mass-acc-ms1 7", "--mass-acc 20 --mass-acc-ms1 7",
                     "--mass-acc 25 --mass-acc-ms1 7", "--mass-acc 23.3 --mass-acc-ms1 7",
                     "--mass-acc 4 --mass-acc-ms1 4"):
            with self.subTest(line), tempfile.TemporaryDirectory() as d:
                out, _ = self._chain(d, _cohort(d))
                with open(os.path.join(out, "window.txt"), "w") as fh:
                    fh.write("7\n")
                with open(os.path.join(out, "massacc.txt"), "w") as fh:
                    fh.write(line + "\n")
                p, argv_log = self._run_step(d, out, "step2_firstpass.sbatch")
                # the fake DIA-NN writes no .quant, so the step still fails afterwards -- what
                # matters is that the guard let it through and DIA-NN got the flags
                self.assertNotIn("outside the plausible band", p.stderr)
                self.assertIn(line, open(argv_log).read())

    def test_after_step1b_the_guard_passes_and_diann_gets_the_pin(self):
        with tempfile.TemporaryDirectory() as d:
            out, _ = self._chain(d, _cohort(d))
            self.assertEqual(self._run_step1b(d, out).returncode, 0)
            p, argv_log = self._run_step(d, out, "step2_firstpass.sbatch")
            self.assertNotIn("step 1b", p.stderr)
            call = open(argv_log).read()
            self.assertIn("--window 7 --mass-acc 20 --mass-acc-ms1 7", call)

    def test_a_failed_mass_accuracy_step1b_stops_and_leaves_no_stale_value(self):
        with tempfile.TemporaryDirectory() as d:
            raws = _cohort(d)
            for r in raws:                                  # no run gets to its MS2 line
                with open(r + ".log", "w") as fh:
                    fh.write(no_ms2(real_log(os.path.basename(r)[:-4])))
            out, _ = self._chain(d, raws)
            for stale, text in (("massacc.txt", "--mass-acc 20 --mass-acc-ms1 7\n"),
                                ("window.txt", "7\n"),
                                ("params.resolved.cfg", "--window 7\n--mass-acc 20\n")):
                with open(os.path.join(out, stale), "w") as fh:   # from an earlier run
                    fh.write(text)
            p = self._run_step1b(d, out)
            log = p.stdout + p.stderr
            self.assertNotEqual(p.returncode, 0, log)
            self.assertNotIn("Traceback", log)
            self.assertIn("FAILED: step 1b measured no scan-window radius and mass accuracy", log)
            self.assertIn("DependencyNeverSatisfied", log)
            for stale in ("massacc.txt", "window.txt", "params.resolved.cfg"):
                self.assertFalse(os.path.exists(os.path.join(out, stale)),
                                 f"a stale {stale} survived a failed step 1b")
            w = json.load(open(os.path.join(out, "window.json")))     # the evidence stays
            self.assertIsNone(w["mass_acc"])

    def test_a_run_without_mass_accuracy_does_not_fail_step1b_when_another_can_replace_it(self):
        with tempfile.TemporaryDirectory() as d:
            raws = _cohort(d) + [_run(d, "Ex01162023_14_TT31", 16 * GB // 10,
                                      log=real_log("Ex01162023_12_TT34")
                                      .replace("Ex01162023_12_TT34", "Ex01162023_14_TT31"))]
            with open(raws[2] + ".log", "w") as fh:          # TT32, the median of these four
                fh.write(no_ms2(real_log("Ex01162023_8_TT32")))
            out, _ = self._chain(d, raws)
            p = self._run_step1b(d, out)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            w = json.load(open(os.path.join(out, "window.json")))
            self.assertEqual(w["failed"], ["Ex01162023_8_TT32.raw"])
            self.assertEqual(w["mass_acc"]["ms2_per_run"], [14.0, 14.0, 17.0])
            self.assertEqual(open(os.path.join(out, "massacc.txt")).read().strip(),
                             "--mass-acc 20 --mass-acc-ms1 7")
            self.assertEqual(w["mass_acc"]["measured_ms2_ppm"], 14.0)


class MassAccPlausibilityTests(unittest.TestCase):
    """A number DIA-NN printed is not automatically a number to pin cohort-wide.

    Until this class existed the only filter anywhere was `> 0`: a cohort whose runs reported MS2
    0.4 / 0.5 / 0.4 ppm pinned `--mass-acc 0.4`, and one reporting 55 / 68 / 60 pinned
    `--mass-acc 60`, both with `"measured": true` and an evidence file beside them. The median
    over three runs does not help -- a wrong FASTA, a wrong species or a miscalibrated batch
    moves all three the same way -- so the pinned value is also checked for magnitude and the
    runs for agreement, and one run is never enough.
    """

    # -- pin_mass_acc: the value that would be pinned -------------------------------------
    def test_a_measured_value_outside_the_orbitrap_band_is_not_pinned(self):
        """MS2 3-30 ppm, MS1 1.5-25 ppm: DIA-NN's own documented tiers span 4-15, it calibrates
        from 25, and its runtime value for a 15k MS2 is 25. Both the too-tight and the too-wide
        cohort below came from running the real probe code."""
        too_tight = probe_window.pin_mass_acc(
            [{"ms2_ppm": 0.4, "ms1_ppm": 7.0}, {"ms2_ppm": 0.5, "ms1_ppm": 7.0},
             {"ms2_ppm": 0.4, "ms1_ppm": 7.0}], documented={"ms1_ppm": 7})
        self.assertIsNone(too_tight["pin_as"], "0.4 ppm reached a DIA-NN command line")
        self.assertTrue(any("outside 3-30 ppm" in r for r in too_tight["rejected"]),
                        too_tight["rejected"])
        too_wide = probe_window.pin_mass_acc(
            [{"ms2_ppm": 55.0, "ms1_ppm": 7.0}, {"ms2_ppm": 68.0, "ms1_ppm": 7.0},
             {"ms2_ppm": 60.0, "ms1_ppm": 7.0}], documented={"ms1_ppm": 7})
        self.assertIsNone(too_wide["pin_as"], "60 ppm reached a DIA-NN command line")
        self.assertTrue(any("60 ppm is outside" in r for r in too_wide["rejected"]),
                        too_wide["rejected"])
        # MS1 has its own band, and it is checked even when MS2 is fine
        ms1_off = probe_window.pin_mass_acc([{"ms2_ppm": 14.0, "ms1_ppm": 0.6},
                                             {"ms2_ppm": 14.0, "ms1_ppm": 0.6},
                                             {"ms2_ppm": 14.0, "ms1_ppm": 0.6}])
        self.assertIsNone(ms1_off["pin_as"])
        self.assertTrue(any(r.startswith("MS1") for r in ms1_off["rejected"]),
                        ms1_off["rejected"])
        # ...and the band does NOT decide the open question: the measured 14 and DIA-NN's own 25
        # for a 15k MS2 are both inside it, and so is the 20 the pilot hand-set.
        for ppm in (14.0, 20.0, 25.0):
            ok = probe_window.pin_mass_acc([{"ms2_ppm": ppm}] * 3, documented={"ms1_ppm": 7})
            self.assertEqual(ok["rejected"], [], f"{ppm} ppm was refused")
            self.assertEqual(ok["measured_ms2_ppm"], ppm)
            # in band and so not refused; the SOP floor then decides which of the two is pinned
            self.assertEqual(ok["pin_as"], f"--mass-acc {max(ppm, 20.0):g} --mass-acc-ms1 7")

    def test_runs_that_disagree_about_a_level_are_not_pinned(self):
        """One acquisition method has one tolerance. A median that no run stands behind is not
        one, however plausible the number looks on its own."""
        split = probe_window.pin_mass_acc([{"ms2_ppm": 14.0, "ms1_ppm": 4.1},
                                           {"ms2_ppm": 14.0, "ms1_ppm": 4.2},
                                           {"ms2_ppm": 30.0, "ms1_ppm": 4.1}])
        self.assertIsNone(split["pin_as"])
        self.assertTrue(any("disagree on MS2" in r for r in split["rejected"]),
                        split["rejected"])
        self.assertAlmostEqual(split["spread"]["ms2_ppm"], 16 / 14.0)
        # the real HIVE cohort, both ways it was measured, still pins
        for per_run in ([14.0, 17.0, 14.0], [14.0, 17.0, 12.0]):
            pin = probe_window.pin_mass_acc([{"ms2_ppm": v} for v in per_run],
                                            documented={"ms1_ppm": 7})
            self.assertEqual(pin["rejected"], [], f"{per_run}: a measured cohort was refused")

    def test_the_median_is_taken_high_so_an_even_count_does_not_pin_the_tighter_value(self):
        """A replaced run makes four probes, and median_low then pins the LOWER of the two middle
        values -- the lossy direction, since a tolerance that is too tight drops identifications
        while one that is too wide only costs specificity. 14 and 20 must pin 20."""
        pin = probe_window.pin_mass_acc([{"ms2_ppm": 14.0}, {"ms2_ppm": 20.0}],
                                        documented={"ms1_ppm": 7})
        self.assertEqual(pin["measured_ms2_ppm"], 20.0)
        self.assertEqual(pin["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")
        # asserted on the MEASURED median throughout: the SOP floor would hide a median_low of
        # 14 behind a pinned 20 and this test would pass without testing anything
        four = probe_window.pin_mass_acc(
            [{"ms2_ppm": 22.0}, {"ms2_ppm": 24.0}, {"ms2_ppm": 26.0}, {"ms2_ppm": 27.0}],
            documented={"ms1_ppm": 7})
        self.assertEqual(four["measured_ms2_ppm"], 26.0)
        self.assertEqual(four["pin_as"], "--mass-acc 26 --mass-acc-ms1 7")
        self.assertIn("median (high)", four["sources"]["ms2_ppm"])
        # odd counts are unchanged: the median is the middle value either way
        odd = probe_window.pin_mass_acc([{"ms2_ppm": 22.0}, {"ms2_ppm": 27.0}, {"ms2_ppm": 22.0}],
                                        documented={"ms1_ppm": 7})
        self.assertEqual(odd["measured_ms2_ppm"], 22.0)

    # -- probe_window.py end to end -------------------------------------------------------
    def _probe(self, d, raws, measure=("window", "mass-acc"), more=(), timeout=60):
        fasta, lib = _fasta_lib(d)
        diann = _exe(os.path.join(d, "diann"), FAKE_DIANN)
        cfg = os.path.join(d, "resolved.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01")
        argv_log = os.path.join(d, "argv.txt")
        argv = [sys.executable, os.path.join(SCRIPTS, "probe_window.py"), "--diann", diann,
                "--raw", *raws, "--fasta", fasta, "--lib", lib, "--threads", "8",
                "--timeout", str(timeout), "--write-cfg", cfg,
                "--workdir", os.path.join(d, "w"), "--measure", *measure, *more]
        p = subprocess.run(argv, capture_output=True, text=True, timeout=240,
                           env=dict(os.environ, DOTNET_ROOT="/opt/fake-dotnet",
                                    FAKE_ARGV_LOG=argv_log))
        return p, cfg, argv_log

    def test_an_implausible_cohort_fails_the_probe_and_writes_nothing(self):
        """Every run agreed, every run logged everything asked, and the answer is still not a
        mass accuracy. The probe must fail exactly as it does for a run that logged nothing: the
        evidence written, the cfg untouched, a non-zero exit, and no `pin_as` anywhere for the
        chain's `json.load(...)['mass_acc']['pin_as']` to pick up."""
        with tempfile.TemporaryDirectory() as d:
            raws = _cohort(d)
            for r in raws:
                with open(r + ".log", "w") as fh:
                    fh.write(ms2_as(real_log(_stem(r)), 60))
            p, cfg, _ = self._probe(d, raws)
            self.assertNotEqual(p.returncode, 0)
            out = json.loads(p.stdout)
            self.assertIsNone(out["mass_acc"])
            self.assertTrue(any("outside 3-30 ppm" in r for r in out["mass_acc_refused"]),
                            out["mass_acc_refused"])
            self.assertIn("Refusing to pin", p.stderr)
            self.assertNotIn("--mass-acc", open(cfg).read())
            self.assertEqual([x["ms2_ppm"] for x in out["probes"]], [60.0, 60.0, 60.0],
                             "the per-run evidence must survive the refusal")
            # the radius travels with it: they are pinned together or not at all
            self.assertIsNone(out["window_radius"])
            self.assertNotIn("--window", open(cfg).read())

    def test_mass_accuracy_is_never_pinned_from_one_run(self):
        """A single surviving run is DIA-NN's own first-run auto mode wearing the evidence of a
        measurement -- and `agree` is True vacuously. Two runs at least, or nothing."""
        with tempfile.TemporaryDirectory() as d:
            raws = _cohort(d)
            for r in raws[1:]:                       # only the smallest run answers
                with open(r + ".log", "w") as fh:
                    fh.write(no_ms2(real_log(_stem(r))))
            p, cfg, _ = self._probe(d, raws)
            self.assertNotEqual(p.returncode, 0)
            out = json.loads(p.stdout)
            good = [x for x in out["probes"] if not x["missing"]]
            self.assertEqual(len(good), 1, [x["missing"] for x in out["probes"]])
            self.assertEqual(good[0]["ms2_ppm"], 14.0, "the one run did measure something")
            self.assertIsNone(out["mass_acc"])
            self.assertTrue(any("at least 2" in r for r in out["mass_acc_refused"]),
                            out["mass_acc_refused"])
            self.assertNotIn("--mass-acc", open(cfg).read())

    def test_one_probe_cannot_be_asked_to_measure_mass_accuracy(self):
        """--max-probes 1 is for the scan window. Asked for mass accuracy it could only ever pin
        from one run, so it is refused before any DIA-NN starts rather than after the wall time."""
        with tempfile.TemporaryDirectory() as d:
            p, cfg, argv_log = self._probe(d, _cohort(d), measure=("mass-acc",),
                                           more=("--max-probes", "1", "--ms1-ppm", "7"))
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("--max-probes", p.stderr)
            self.assertFalse(os.path.exists(argv_log), "DIA-NN ran for a pin it could not make")
            self.assertNotIn("--mass-acc", open(cfg).read())
            # the window alone is still a one-probe job
            p, _, _ = self._probe(d, _cohort(d), measure=("window",),
                                  more=("--max-probes", "1"))
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual(json.loads(p.stdout)["window_radius"], 7)

    def test_a_documented_level_outside_the_band_is_refused_before_diann_starts(self):
        """--ms1-ppm/--ms2-ppm are pinned as given and never measured, so pin_mass_acc cannot
        check them. They are checked as they are parsed, so every number that can reach `pin_as`
        has been through the band."""
        with tempfile.TemporaryDirectory() as d:
            for flag, value in (("--ms1-ppm", "999"), ("--ms1-ppm", "0.2"),
                                ("--ms2-ppm", "120")):
                with self.subTest(flag=flag, value=value):
                    p, cfg, argv_log = self._probe(d, _cohort(d), measure=("mass-acc",),
                                                   more=(flag, value))
                    self.assertNotEqual(p.returncode, 0)
                    self.assertIn("outside", p.stderr)
                    self.assertIn(flag, p.stderr)
                    self.assertFalse(os.path.exists(argv_log), "DIA-NN ran on a junk tier")
            p, _, _ = self._probe(d, _cohort(d), measure=("mass-acc",),
                                  more=("--ms1-ppm", "7"))
            self.assertEqual(p.returncode, 0, p.stderr)


class SopFloorTests(unittest.TestCase):
    """A measured level is pinned at max(measured, SOP): the measurement is used only where it is
    WIDER than the facility's SOP tolerance.

    Where the probe earns its keep is an instrument that genuinely needs a wider window than the
    SOP -- nothing else would catch that. A measured value TIGHTER than the SOP buys nothing and
    costs identifications: on the one cohort benchmarked (references/diann_parallel.md) the
    measured 14/7 gave 18,476 precursors against 19,592 at the SOP's 20/7, same runs, same
    library, same FDR. Flooring keeps the win and drops the loss.

    The floor must never erase the measurement -- that would be the same false provenance the
    band exists to prevent -- and must never rescue a measurement the band refused.
    """

    _probe = MassAccPlausibilityTests._probe        # the same probe_window.py CLI runner

    def test_the_sop_floor_has_exactly_one_definition(self):
        """`probe_window` does not re-type 20/7: it imports the SOP from the module that already
        owns every mass-accuracy table, so the floor moves when the SOP does. Nothing else in the
        skill defines an SOP tolerance -- make_presets.py's 20/20 is Radiant's and FragPipe's
        VENDOR default, which that file itself calls too wide for narrow-window data."""
        import estimate_params
        self.assertIs(probe_window.SOP_MASS_ACC, estimate_params.SOP_MASS_ACC)
        import diann_parallel
        self.assertEqual(diann_parallel.SOP_MASS_ACC_FLAGS,
                         {"--mass-acc": estimate_params.SOP_MASS_ACC["ms2_ppm"],
                          "--mass-acc-ms1": estimate_params.SOP_MASS_ACC["ms1_ppm"]})
        # the floor tracks the SOP rather than a literal
        src = open(os.path.join(SCRIPTS, "probe_window.py")).read()
        floor = src[src.index("        floor = float("):src.index("        floor = float(") + 120]
        self.assertIn("SOP_MASS_ACC[key]", floor, floor)

    def test_a_measurement_wider_than_the_sop_is_the_one_that_is_pinned(self):
        """The case the probe exists for: this instrument needs more than the SOP allows."""
        pin = probe_window.pin_mass_acc([{"ms2_ppm": 25.0}, {"ms2_ppm": 26.0},
                                         {"ms2_ppm": 25.0}], documented={"ms1_ppm": 7})
        self.assertEqual(pin["measured_ms2_ppm"], 25.0)
        self.assertEqual(pin["pinned_ms2_ppm"], 25.0)
        self.assertFalse(pin["floored"]["ms2_ppm"])
        self.assertEqual(pin["pin_as"], "--mass-acc 25 --mass-acc-ms1 7")
        self.assertIn("at or above the 20 ppm SOP floor", pin["sources"]["ms2_ppm"])
        # exactly at the floor is not floored either -- max() of equals
        at = probe_window.pin_mass_acc([{"ms2_ppm": 20.0}] * 3, documented={"ms1_ppm": 7})
        self.assertFalse(at["floored"]["ms2_ppm"])
        self.assertEqual(at["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")

    def test_a_measurement_tighter_than_the_sop_pins_the_sop_and_records_both(self):
        """The validation cohort's own case: measured 14, searched at 20. A reader of the
        evidence must be able to see BOTH numbers and that the floor is why they differ."""
        pin = probe_window.pin_mass_acc([{"ms2_ppm": 14.0, "ms1_ppm": 4.1},
                                         {"ms2_ppm": 17.0, "ms1_ppm": 4.3},
                                         {"ms2_ppm": 14.0, "ms1_ppm": 4.2}])
        self.assertEqual(pin["measured_ms2_ppm"], 14.0)
        self.assertEqual(pin["pinned_ms2_ppm"], 20.0)
        self.assertEqual(pin["measured_ms1_ppm"], 4.2)
        self.assertEqual(pin["pinned_ms1_ppm"], 7.0)
        self.assertEqual(pin["floored"], {"ms2_ppm": True, "ms1_ppm": True})
        self.assertEqual(pin["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")
        self.assertEqual(pin["ms2_per_run"], [14.0, 17.0, 14.0], "the per-run values survive")
        self.assertEqual(pin["sop_floor"], {"ms1_ppm": 7.0, "ms2_ppm": 20.0})
        for key in ("ms2_ppm", "ms1_ppm"):
            self.assertIn("TIGHTER than the SOP floor", pin["sources"][key])
        self.assertIn("measured 14 ppm", pin["sources"]["ms2_ppm"])
        self.assertIn("pinned at 20 ppm", pin["sources"]["ms2_ppm"])

    def test_the_floor_does_not_rescue_a_measurement_the_band_refused(self):
        """The band and the floor answer different questions. 0.4 ppm is tighter than the SOP,
        but it is not a mass accuracy at all -- flooring it to 20 would turn the probe's one
        unambiguous failure signal into a silent fallback to the SOP, with `measured: true` and
        an evidence file to back it up. It stays a probe FAILURE, and `floored` stays false."""
        for per_run in ([0.4, 0.5, 0.4], [1.0, 1.0, 1.0]):
            with self.subTest(per_run):
                pin = probe_window.pin_mass_acc([{"ms2_ppm": v} for v in per_run],
                                                documented={"ms1_ppm": 7})
                self.assertIsNone(pin["pin_as"])
                self.assertTrue(pin["rejected"])
                self.assertFalse(pin["floored"]["ms2_ppm"], "an implausible value was floored up")
                self.assertEqual(pin["pinned_ms2_ppm"], pin["measured_ms2_ppm"],
                                 "the SOP became a fallback for a refused measurement")
        # ...and the same for a level refused on spread rather than magnitude
        split = probe_window.pin_mass_acc([{"ms2_ppm": 4.0}, {"ms2_ppm": 4.0},
                                           {"ms2_ppm": 14.0}], documented={"ms1_ppm": 7})
        self.assertIsNone(split["pin_as"])
        self.assertFalse(split["floored"]["ms2_ppm"])

    def test_the_probe_writes_the_floored_value_and_says_both_numbers(self):
        """End to end: the cfg and the job log, not only the JSON."""
        with tempfile.TemporaryDirectory() as d:
            p, cfg, _ = self._probe(d, _cohort(d), measure=("mass-acc",),
                                    more=("--ms1-ppm", "7"))
            self.assertEqual(p.returncode, 0, p.stderr)
            ma = json.loads(p.stdout)["mass_acc"]
            self.assertEqual((ma["measured_ms2_ppm"], ma["pinned_ms2_ppm"]), (14.0, 20.0))
            self.assertIn("--mass-acc 20\n", open(cfg).read())
            self.assertNotIn("--mass-acc 14\n", open(cfg).read())
            self.assertIn("MS2: measured 14 ppm, PINNED 20 ppm -- the SOP floor", p.stderr)

    def test_a_cohort_that_needs_a_wider_window_gets_it_end_to_end(self):
        """The floor must not flatten every cohort to the SOP. A cohort whose runs really do
        measure wider than the SOP searches at what it measured."""
        with tempfile.TemporaryDirectory() as d:
            raws = _cohort(d)
            for r in raws:
                with open(r + ".log", "w") as fh:
                    fh.write(ms2_as(real_log(_stem(r)), 26))
            p, cfg, _ = self._probe(d, raws, measure=("mass-acc",), more=("--ms1-ppm", "7"))
            self.assertEqual(p.returncode, 0, p.stderr)
            ma = json.loads(p.stdout)["mass_acc"]
            self.assertEqual((ma["measured_ms2_ppm"], ma["pinned_ms2_ppm"]), (26.0, 26.0))
            self.assertFalse(ma["floored"]["ms2_ppm"])
            self.assertIn("--mass-acc 26\n", open(cfg).read())
            self.assertNotIn("PINNED", p.stderr, "a floor was announced where none applied")

def _stem(raw):
    """'/tmp/x/Ex01162023_12_TT34.raw' -> 'Ex01162023_12_TT34'."""
    return os.path.basename(raw)[:-len(".raw")]


class BatchDriftGuidanceTests(unittest.TestCase):
    """What step 1b pins is the instrument's calibration on the days the probed runs were
    acquired. A cohort acquired over a week can cross a recalibration, a cleaning or a lock-mass
    change, and three representative runs cannot see it -- the spread check only catches drift
    big enough to move the runs it happened to probe. The reference had no guidance at all, so a
    user with a multi-day batch had nothing telling them the single pinned value assumes one
    measurement condition."""

    REF = os.path.join(os.path.dirname(HERE), "references", "diann_parallel.md")

    def test_the_reference_records_the_floor_decision_not_an_open_question(self):
        """The benchmark table's "Open question" was the maintainer's to settle, and it was:
        measure, but floor at the SOP. The reference has to record the decision and why, or the
        next reader re-opens it. The measured numbers in the table stay as they are."""
        text = " ".join(open(self.REF).read().split())
        self.assertNotIn("Open question", text, "the question is settled; the doc still asks it")
        self.assertIn("max(measured, SOP)", text)
        self.assertIn("SOP_MASS_ACC", text)
        for kept in ("18,476", "19,592", "2,618", "25 ppm"):
            self.assertIn(kept, text, f"the benchmark lost {kept}")
        where = text.find("Settled \u2014 measure")       # the heading, not the cross-reference
        self.assertNotEqual(where, -1, "the decision has no section of its own")
        para = text[where:where + 2200]
        self.assertIn("wider", para, "no reason given for flooring rather than pinning")
        self.assertIn("measured_ms2_ppm", para, "the doc does not say where both numbers are")
        self.assertIn("never floored", para, "a documented level's exemption is not recorded")
    def test_the_reference_says_to_split_a_long_batch(self):
        text = " ".join(open(self.REF).read().lower().split())   # the doc hard-wraps at 100
        where = text.find("mass accuracy drifts")
        self.assertNotEqual(where, -1, "no guidance on mass-accuracy drift within a batch")
        para = text[where:where + 1600]
        self.assertIn("split", para)
        self.assertIn("window.json", para, "no way given to compare the parts")
        for word in ("recalibration", "acquisition date"):
            self.assertIn(word, para, word)

if __name__ == "__main__":
    unittest.main(verbosity=2)
