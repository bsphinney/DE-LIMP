#!/usr/bin/env python3
"""
The pre-search probe must survive a stale NFS handle, and a probe that measures nothing must not
take the search down with it.

fran-5b, 2026-09-29 17:00 to 09-30 08:00: 62 of 261 step-1b jobs (24%), on 28 of 56 nodes, died in
probe_window.run_probe with `OSError: [Errno 116] Stale file handle` at `pending += src.read()`.
The probe tailed <workdir>/probe.log -- on Flinders NFS -- through one handle held open for the
whole ~20-minute DIA-NN run. DIA-NN itself succeeded; afterok then left steps 2-5 (183 jobs)
DependencyNeverSatisfied.

  1. The live log is on node-local storage; <workdir>/probe.log is a copy, published on every
     way out of run_probe.
  2. The tail reopens on ESTALE at its own byte offset (_Tail); an error that still gets through
     fails that RUN (io_error, replaced like any run that logged nothing), not the probe.
  3. Only when the probe's OWN MACHINERY failed -- its log unreadable, a crash, a time limit --
     does the job go on without a measurement: it retries a crash or an unreadable log once,
     then falls back (probe_fallback.py): window `auto` (DIA-NN chooses per run) and mass
     accuracy at the documented level / the SOP tagged DEFAULT, exiting 0 so steps 2-5 run --
     recorded as a fallback everywhere, with Methods that say the measurement failed (DE-LIMP
     rule 2) and a CAUTION wherever the results are read. Everything else still fails the job,
     never retried or fallen back past (dda-review, 2026-09-30): a refused mass accuracy (any
     run's logged value, complete or not), the environment (no .NET), the probe's arguments,
     runs DIA-NN finished without logging it, a signal.

The chain tests inject the ESTALE into the probe's real process (estale_inject.py): any read of
a file named probe.log, wherever it is -- so the same tests crash the probe of skill <= 979dfd5.
"""
import errno
import glob
import json
import os
import stat
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
from job_env import job_env  # noqa: E402  (env for running job scripts)

import probe_window  # noqa: E402
import probe_fallback  # noqa: E402
import diann_parallel as dp  # noqa: E402
import make_methods  # noqa: E402
import make_analysis_html  # noqa: E402
import record_run  # noqa: E402
from estimate_params import SOP_DEFAULT_PUBLISHED  # noqa: E402
import test_step1b_mass_accuracy_probe as fx  # noqa: E402  (captured logs and the probe fake)
from estale_inject import estale_env  # noqa: E402

MB = 1024 ** 2
STALE = os.strerror(errno.ESTALE)

# Logs its argv; logs a radius at once, then "searches" for FAKE_SLEEP seconds. FAKE_NO_DOTNET:
# the line DIA-NN prints, and exits 0 on, when it cannot read .raw without .NET 8.
FAKE_DIANN = r"""#!/bin/bash
[ -n "${FAKE_ARGV_LOG:-}" ] && echo "$@" >> "$FAKE_ARGV_LOG"
f=""
while [ $# -gt 0 ]; do case "$1" in --f) f="$2"; shift;; esac; shift; done
# the data answering: DIA-NN finishes the run without logging a radius
[ -n "${FAKE_NO_RADIUS:-}" ] && { echo "fake DIA-NN: $f"; echo "Finished"; exit 0; }
# goes silent before the radius: the probe's --timeout
[ -n "${FAKE_HANG:-}" ] && { echo "fake DIA-NN: $f"; exec sleep 30; }
if [ -n "${FAKE_NO_DOTNET:-}" ]; then
  echo "ERROR: cannot read .raw files, please download and install .NET Runtime 8: 8.0.17 or later https://dotnet.microsoft.com/en-us/download/dotnet/8.0 : 1"
  exit 0
fi
echo "fake DIA-NN: $f"
sleep "${FAKE_FIRST:-0}"
echo "[0:01] Loading run $f"
sleep "${FAKE_FIRST:-0}"
echo "[0:04] Scan window radius set to 7"
exec sleep "${FAKE_SLEEP:-20}"
"""

def _exe(path, body):
    with open(path, "w") as fh:
        fh.write(body)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def _read(path):
    with open(path) as fh:
        return fh.read()


class _Flaky:
    """A real handle whose read() raises `err` on the reads numbered in `fail_on` (counted over
    every handle of the file, from 1) WITHOUT moving the file position -- a failed NFS read
    need not have."""

    def __init__(self, fh, state):
        self.fh, self.state = fh, state

    def read(self, *a):
        self.state["reads"] += 1
        if self.state["reads"] in self.state["fail_on"]:
            raise self.state["err"]
        return self.fh.read(*a)

    def seek(self, *a):
        return self.fh.seek(*a)

    def close(self):
        return self.fh.close()


def flaky_open(fail_on, err=None, match="probe_log_"):
    """(fake open for probe_window, its state): reads of files whose path contains `match`,
    opened "rb", fail as _Flaky says. Every such path is recorded in state["opened"]."""
    state = {"reads": 0, "fail_on": set(fail_on), "opened": [],
             "err": err or OSError(errno.ESTALE, STALE)}

    def _open(path, mode="r", *a, **k):
        fh = open(path, mode, *a, **k)
        if mode == "rb" and match in str(path):
            state["opened"].append(str(path))
            return _Flaky(fh, state)
        return fh
    return _open, state


class TailTests(unittest.TestCase):
    """_Tail: reopen on ESTALE at its OWN offset -- once."""

    def test_estale_mid_tail_reopens_at_its_own_offset(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "probe_log_x.log")
            with open(path, "wb") as fh:
                fh.write(b"line1\n")
            fake, st = flaky_open({2})
            with mock.patch.object(probe_window, "open", fake, create=True):
                t = probe_window._Tail(path)
                got = t.read()
                with open(path, "ab") as fh:
                    fh.write(b"line2\n")
                got += t.read()                 # read 2 raises ESTALE; read 3 is the reopen
                with open(path, "ab") as fh:
                    fh.write(b"line3\n")
                got += t.read()
                t.close()
            self.assertEqual(got, b"line1\nline2\nline3\n", "lost or repeated a line")
            self.assertEqual(len(st["opened"]), 2, "the stale handle was not reopened")
            self.assertEqual(t.offset, len(got))

    def test_a_second_estale_in_a_row_is_raised(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "probe_log_x.log")
            with open(path, "wb") as fh:
                fh.write(b"line1\n")
            fake, st = flaky_open({1, 2})
            with mock.patch.object(probe_window, "open", fake, create=True):
                t = probe_window._Tail(path)
                with self.assertRaises(OSError) as cm:
                    t.read()
            self.assertEqual(cm.exception.errno, errno.ESTALE)
            self.assertEqual(len(st["opened"]), 2)
            self.assertIsNone(t.fh)

    def test_other_errors_are_not_retried(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "probe_log_x.log")
            with open(path, "wb") as fh:
                fh.write(b"line1\n")
            fake, st = flaky_open({1}, err=OSError(errno.EIO, os.strerror(errno.EIO)))
            with mock.patch.object(probe_window, "open", fake, create=True):
                with self.assertRaises(OSError):
                    probe_window._Tail(path).read()
            self.assertEqual(len(st["opened"]), 1)


class RunProbeLogTests(unittest.TestCase):
    """run_probe against a fake DIA-NN: the live log is node-local, the workdir gets a copy."""

    def _probe(self, d, timeout=60, diann=None, env=None):
        diann = diann or _exe(os.path.join(d, "diann"), FAKE_DIANN)
        wd = os.path.join(d, "nfs_workdir", "probe1")
        local = os.path.join(d, "node_local")
        os.makedirs(local, exist_ok=True)
        with mock.patch.object(probe_window.tempfile, "tempdir", local), \
                mock.patch.dict(os.environ, dict({"FAKE_FIRST": "0.3"}, **(env or {}))):
            t0 = time.time()
            r = probe_window.run_probe(diann, os.path.join(d, "run1.mzML"),
                                       os.path.join(d, "db.fasta"), os.path.join(d, "lib"),
                                       4, timeout, workdir=wd, measure=("window",))
        return r, wd, local, time.time() - t0

    def test_estale_mid_tail_still_reads_the_radius_and_publishes_the_log(self):
        with tempfile.TemporaryDirectory() as d:
            fake, st = flaky_open({2})
            with mock.patch.object(probe_window, "open", fake, create=True):
                r, wd, local, _ = self._probe(d)
            self.assertEqual(r["radius"], 7, r["lines"])
            self.assertIsNone(r["io_error"])
            self.assertIsNone(r["log_publish_error"])
            self.assertGreaterEqual(len(st["opened"]), 2, "the ESTALE was never hit")
            for p in st["opened"]:                  # tailed on node-local storage, not the workdir
                self.assertTrue(p.startswith(local + os.sep), p)
            self.assertEqual(r["log"], os.path.join(wd, "probe.log"))
            log = _read(r["log"])
            self.assertEqual(log.count("Scan window radius set to 7"), 1, log)
            self.assertEqual(log.count("Loading run"), 1, log)
            self.assertEqual(os.listdir(local), [], "the node-local log was left behind")

    def test_a_log_that_stays_unreadable_fails_the_run_not_the_probe(self):
        with tempfile.TemporaryDirectory() as d:
            fake, _ = flaky_open(range(1, 1000))
            with mock.patch.object(probe_window, "open", fake, create=True):
                r, wd, local, secs = self._probe(d)
            self.assertIsNone(r["radius"])
            self.assertIn(STALE, r["io_error"])
            self.assertEqual(r["missing"], ["radius"])
            self.assertFalse(r["environmental"], "one run's I/O error must not stop the probe")
            self.assertLess(secs, 15, "DIA-NN was left running after the log went unreadable")
            self.assertIn("could not be read", _read(os.path.join(wd, "probe.log")))
            self.assertEqual(os.listdir(local), [])

    def test_every_way_out_of_run_probe_publishes_the_log(self):
        with tempfile.TemporaryDirectory() as d:
            r, wd, _, _ = self._probe(d, timeout=0)                 # no time left to start
            self.assertIn("no time left", _read(os.path.join(wd, "probe.log")))
        with tempfile.TemporaryDirectory() as d:
            r, wd, _, _ = self._probe(d, diann=os.path.join(d, "no-such-diann"))
            self.assertTrue(r["environmental"])
            self.assertIn("could not start DIA-NN", _read(os.path.join(wd, "probe.log")))
        with tempfile.TemporaryDirectory() as d:
            r, wd, _, _ = self._probe(d)
            self.assertEqual(r["radius"], 7)
            self.assertIn("Scan window radius set to 7", _read(os.path.join(wd, "probe.log")))

    def test_a_log_that_cannot_be_published_does_not_lose_the_measurement(self):
        with tempfile.TemporaryDirectory() as d:
            with mock.patch.object(probe_window.shutil, "copyfile",
                                   side_effect=OSError(errno.ESTALE, STALE)):
                r, wd, local, _ = self._probe(d)
            self.assertEqual(r["radius"], 7)
            self.assertIn(STALE, r["log_publish_error"])
            self.assertEqual(os.listdir(local), [])


class ProbeAttemptsBashTests(unittest.TestCase):
    """diann_parallel.probe_attempts(), run in bash under `set -euo pipefail` (as a job runs it)
    against a probe that exits as told: RC_<n> for attempt n."""

    FAKE_PROBE = r"""#!/bin/bash
n=$(( $(cat "$COUNT" 2>/dev/null || echo 0) + 1 )); echo "$n" > "$COUNT"
echo "$*" >> "$ARGS"
rc_var="RC_$n"; rc=${!rc_var:-0}
if [ "$rc" -eq 0 ]; then echo '{"window_radius": 7, "stopped_because": "measured"}'
else echo '{"window_radius": null, "stopped_because": "max_failures", "failed": ["a.mzML"],
  "probes": [{"io_error": "OSError: [Errno 116] Stale file handle"}]}'; fi
exit "$rc"
"""

    def _run(self, d, rcs, budget=None):
        fake = _exe(os.path.join(d, "fake_probe"), self.FAKE_PROBE)
        ev, wtxt, fb = (os.path.join(d, n) for n in ("window.json", "window.txt", "fb.json"))
        kw = dict(reset=["echo reset >> \"$RESETS\""], window_file=wtxt, fallback_out=fb)
        if budget is None:
            lines = dp.probe_attempts(f"bash {fake}", "--qvalue 0.01", ev, ["window"], {}, **kw)
        else:
            with mock.patch.object(dp, "PROBE_BUDGET_S", budget):
                lines = dp.probe_attempts(f"bash {fake}", "--qvalue 0.01", ev, ["window"], {},
                                          **kw)
        script = os.path.join(d, "attempts.sh")
        with open(script, "w") as fh:
            fh.write("set -euo pipefail\n" + "\n".join(lines)
                     + '\necho "RC=$PROBE_RC FB=$PROBE_FALLBACK N=$PROBE_ATTEMPTS"\n')
        env = job_env(d, COUNT=os.path.join(d, "count"), ARGS=os.path.join(d, "args"),
                      RESETS=os.path.join(d, "resets"),
                      **{f"RC_{i}": rc for i, rc in enumerate(rcs, 1)})
        p = subprocess.run(["bash", script], capture_output=True, text=True, env=env,
                           timeout=60)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        state = p.stdout.strip().splitlines()[-1]
        return p, state, {"window.txt": wtxt, "fb": fb, "ev": ev,
                          "args": os.path.join(d, "args"), "resets": os.path.join(d, "resets")}

    def test_measured_first_time(self):
        with tempfile.TemporaryDirectory() as d:
            p, state, f = self._run(d, [0])
            self.assertEqual(state, "RC=0 FB=0 N=1")
            self.assertFalse(os.path.exists(f["window.txt"]))
            self.assertFalse(os.path.exists(f["resets"]))

    def test_a_crash_or_an_unreadable_log_is_retried_and_a_retry_that_measures_needs_nothing(self):
        for rc in (probe_window.EXIT_CRASH, probe_window.EXIT_IO_ERROR):
            with self.subTest(rc=rc), tempfile.TemporaryDirectory() as d:
                p, state, f = self._run(d, [rc, 0])
                self.assertEqual(state, "RC=0 FB=0 N=2")
                self.assertIn("retrying once", p.stderr)
                self.assertEqual(_read(f["resets"]).split(), ["reset"])
                self.assertTrue(os.path.exists(f["ev"] + ".attempt1"),
                                "attempt 1's evidence lost")
                self.assertFalse(os.path.exists(f["fb"]))
                # ONE budget, set before both attempts: the retry gets what is left of it
                b1, b2 = (int(ln.split("--budget ")[1].split()[0])
                          for ln in _read(f["args"]).splitlines())
                self.assertLessEqual(b2, b1)
                self.assertLessEqual(b1, dp.PROBE_BUDGET_S)

    def test_two_machinery_failures_fall_back(self):
        for rcs in ([1, 1], [7, 7], [1, 7]):
            with self.subTest(rcs=rcs), tempfile.TemporaryDirectory() as d:
                p, state, f = self._run(d, rcs)
                self.assertEqual(state, f"RC={rcs[1]} FB=1 N=2")
                self.assertEqual(_read(f["window.txt"]).strip(), "auto")
                rec = json.loads(_read(f["fb"]))
                self.assertIn("I/O error on the probe log", rec["reason"])
                self.assertIn(f"2 attempts, last exit {rcs[1]}", rec["reason"])
                self.assertIn("FALLBACK", p.stderr)

    def test_a_time_limit_falls_back_without_a_retry(self):
        with tempfile.TemporaryDirectory() as d:
            p, state, f = self._run(d, [probe_window.EXIT_TIMED_OUT, 0])
            self.assertEqual(state, "RC=8 FB=1 N=1")
            self.assertNotIn("retrying", p.stderr)
            self.assertEqual(_read(f["window.txt"]).strip(), "auto")

    def test_nothing_else_is_ever_retried_or_fallen_back_past(self):
        """dda-review, 2026-09-30: a refusal, the environment (no .NET), argparse and the
        probe's own argument checks, runs DIA-NN finished without logging it, a signal. No
        fallback can fix them, and one would record a generator bug or a wrong FASTA as a
        harmless "fallback"."""
        for rc in (probe_window.EXIT_USAGE, probe_window.EXIT_REFUSED,
                   probe_window.EXIT_ENVIRONMENT, probe_window.EXIT_CONFIG,
                   probe_window.EXIT_NOT_MEASURED, 143, 127):
            with self.subTest(rc=rc), tempfile.TemporaryDirectory() as d:
                p, state, f = self._run(d, [rc, 0])
                self.assertEqual(state, f"RC={rc} FB=0 N=1")
                self.assertNotIn("retrying", p.stderr)
                self.assertFalse(os.path.exists(f["window.txt"]))
                self.assertFalse(os.path.exists(f["fb"]))

    def test_the_exit_statuses_are_the_rule(self):
        self.assertEqual(probe_window.RETRY_ON, (1, 7))
        self.assertEqual(probe_window.FALLBACK_ON, (1, 7, 8))
        # every deliberate stop explains itself in a failed job's log
        self.assertEqual(set(probe_window.EXIT_MEANING), {2, 3, 4, 5, 6})

    def test_no_retry_without_the_time_for_one(self):
        with tempfile.TemporaryDirectory() as d:
            p, state, f = self._run(d, [1, 0], budget=100)
            self.assertEqual(state, "RC=1 FB=1 N=1")
            self.assertIn("not retried", p.stderr)
            self.assertIn("1 attempt,", json.loads(_read(f["fb"]))["reason"])


class FallbackRecordTests(unittest.TestCase):
    """probe_fallback.py's record, and what the Methods and the run record make of it."""

    def test_window_and_mass_accuracy_fall_back_to_auto_documented_and_the_sop(self):
        rec = probe_fallback.fallback(["window", "mass-acc"], {"ms1_ppm": 7.0}, "why")
        self.assertEqual(rec["window"]["value"], "auto")
        self.assertTrue(rec["window"]["source"].startswith("fallback (probe failed: why)"))
        ma = rec["mass_acc"]
        self.assertEqual(ma["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")
        self.assertFalse(ma["measured"])
        self.assertEqual(ma["documented"], {"--mass-acc-ms1": 7.0})
        self.assertEqual(ma["default"], {"--mass-acc": 20.0})
        self.assertIn("DEFAULT, not user-confirmed", ma["source"])
        # nothing documented: both levels at the SOP, both DEFAULT
        ma = probe_fallback.fallback(["mass-acc"], {}, "why")["mass_acc"]
        self.assertEqual(ma["default"], {"--mass-acc": 20.0, "--mass-acc-ms1": 7.0})

    def test_a_probe_that_left_no_evidence_is_said_to_have_crashed(self):
        with tempfile.TemporaryDirectory() as d:
            why = probe_fallback.reason_from(os.path.join(d, "missing.json"), 1, 2)
            self.assertIn("without writing its evidence", why)
            self.assertIn("2 attempts", why)

    def test_the_script_writes_every_file_and_records_the_fallback_in_provenance(self):
        with tempfile.TemporaryDirectory() as d:
            prov = os.path.join(d, "search_provenance.json")
            with open(prov, "w") as fh:
                json.dump({"engine": "diann", "scan_window": {"source": "measured at run time"},
                           "result": {"scan_window": {"source": "measured at run time"},
                                      "mass_acc": {"measured": True}}}, fh)
            cfg = os.path.join(d, "params.resolved.cfg.tmp")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01")                 # no trailing newline
            files = {n: os.path.join(d, n) for n in ("window.txt", "massacc.txt", "fb.json")}
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "probe_fallback.py"),
                                "--measure", "window", "mass-acc", "--ms1-ppm", "7",
                                "--exit-code", "1", "--attempts", "2",
                                "--evidence", os.path.join(d, "window.json"),
                                "--window-file", files["window.txt"],
                                "--massacc-file", files["massacc.txt"], "--write-cfg", cfg,
                                "--provenance", prov, "--out", files["fb.json"]],
                               capture_output=True, text=True, timeout=60)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual(_read(files["window.txt"]), "auto\n")
            self.assertEqual(_read(files["massacc.txt"]), "--mass-acc 20 --mass-acc-ms1 7\n")
            self.assertEqual(dp.cfg_tokens(cfg), ["--qvalue", "0.01", "--mass-acc", "20",
                                                  "--mass-acc-ms1", "7"])
            pv = json.loads(_read(prov))
            self.assertTrue(pv["probe_fallback"]["fallback"])
            self.assertTrue(pv["scan_window"]["fallback"])
            self.assertIsNone(pv["scan_window"]["value"])
            self.assertEqual(pv["result"]["scan_window"], pv["scan_window"])
            self.assertFalse(pv["result"]["mass_acc"]["measured"])
            self.assertEqual(pv["result"]["mass_acc"]["value_file"], files["massacc.txt"])
            self.assertIn("WARNING: the probe measured nothing", p.stderr)

    def _search_record(self, d, measure):
        cfg = os.path.join(d, "params.resolved.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01\n--cut K*,R*\n"
                     + ("--mass-acc 20\n--mass-acc-ms1 7\n" if "mass-acc" in measure else ""))
        fb = probe_fallback.fallback(measure, {"ms1_ppm": 7.0}, "the runs tried did not log it")
        prov = {"engine": "diann", "version": "2.7.0", "resolved_params_file": cfg,
                "search_mode": "parallel_5step", "probe_fallback": fb, "result": {}}
        if "mass_acc" in fb:
            prov["result"]["mass_acc"] = fb["mass_acc"]
        path = os.path.join(d, "search_provenance.json")
        with open(path, "w") as fh:
            json.dump(prov, fh)
        return make_methods.search_record(search_prov=path)

    def test_the_methods_say_the_window_was_set_by_dia_nn_because_the_measurement_failed(self):
        with tempfile.TemporaryDirectory() as d:
            rec = self._search_record(d, ["window"])
            para = make_methods.search_paragraph(rec)
            self.assertIn("The scan window was to be measured on representative runs before the "
                          "search, but that measurement failed; DIA-NN set the scan window "
                          "automatically, for each run.", para)
            self.assertNotIn("mass tolerances given here were not measured", para)

    def test_the_methods_say_a_fallen_back_mass_accuracy_was_not_measured(self):
        with tempfile.TemporaryDirectory() as d:
            rec = self._search_record(d, ["window", "mass-acc"])
            para = make_methods.search_paragraph(rec)
            self.assertIn("The scan window and the mass accuracy were to be measured", para)
            self.assertIn("DIA-NN set the scan window automatically, for each run, and the mass "
                          "tolerances given here were not measured on these data.", para)
            # the SOP level carries the DEFAULT tag; the documented one does not
            self.assertIn(f"fragment (MS2) 20 ppm ({SOP_DEFAULT_PUBLISHED})", para)
            self.assertIn("precursor (MS1) 7 ppm and", para)
            self.assertEqual(rec["ms2_tol"]["default"], SOP_DEFAULT_PUBLISHED)
            self.assertIsNone(rec["ms1_tol"].get("default"))
            # the run record's table says FALLBACK, never "measured"
            table = "\n".join(record_run.render_parameters({"parameters": {
                "record": rec, "why": {}, "mass_accuracy": {
                    "massacc_txt": "--mass-acc 20 --mass-acc-ms1 7",
                    "measured": rec and {"measured": False}},
                "scan_window": {"value": None, "fallback": True, "source": "fallback (...)"}}}))
            self.assertIn("auto -- DIA-NN chose it per run (FALLBACK, not measured)", table)
            self.assertIn("massacc.txt -- FALLBACK, NOT measured", table)
            self.assertNotIn("Mass accuracy measured", table)

    def test_no_fallback_no_sentence(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n--window 7\n--mass-acc 20\n--mass-acc-ms1 7\n")
            path = os.path.join(d, "search_provenance.json")
            with open(path, "w") as fh:
                json.dump({"engine": "diann", "version": "2.7.0", "resolved_params_file": cfg,
                           "result": {}}, fh)
            rec = make_methods.search_record(search_prov=path)
            self.assertIsNone(rec["probe_fallback"])
            self.assertNotIn("measurement failed", make_methods.search_paragraph(rec))


class Step1bEstaleChainTests(unittest.TestCase):
    """The generated step1b_window.sbatch, run with bash as a compute node runs it, with the
    ESTALE injected into the probe's own process (estale_inject.py): any read of a probe.log,
    wherever it is. With skill <= 979dfd5 in place of the fix, the transient case crashes the
    probe with the incident's traceback and fails step 1b."""

    def _chain(self, d, n=5):
        runs = []
        for i in range(n):
            p = os.path.join(d, f"run{i}.mzML")
            with open(p, "wb") as fh:
                fh.truncate((30 + i) * MB)
            runs.append(p)
        fasta = os.path.join(d, "db.fasta")
        with open(fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDER\n")
        cfg = os.path.join(d, "p.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01\n--mass-acc 20\n--mass-acc-ms1 7\n")
        diann = _exe(os.path.join(d, "diann"), FAKE_DIANN)
        out = os.path.join(d, "out")
        g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                            "--diann", diann, "--raw", *runs, "--fasta", fasta, "--out", out,
                            "--cfg", cfg, "--threads-per-file", "4"],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(g.returncode, 0, g.stderr)
        info = json.loads(g.stdout)
        with open(os.path.join(out, "step1.predicted.speclib"), "w") as fh:
            fh.write("lib")                             # stands in for step 1's library
        with open(os.path.join(out, "search_provenance.json"), "w") as fh:
            # as run_search.py writes it when it generates the chain
            json.dump({"engine": "diann", "scan_window": info["scan_window"],
                       "mass_acc": info["mass_acc"], "result": info}, fh)
        return out

    def _run(self, d, out, name, estale=None, **extra):
        env = job_env(d, base={k: v for k, v in os.environ.items() if k != "DOTNET_ROOT"},
                      **(estale_env(d, estale) if estale else {}), **extra)
        return subprocess.run(["bash", os.path.join(out, name)], cwd=out, capture_output=True,
                              text=True, env=env, timeout=240)

    def _failed_without_fallback(self, out, p, status):
        log = p.stdout + p.stderr
        self.assertNotEqual(p.returncode, 0, log)
        self.assertNotIn("Traceback", log)
        self.assertNotIn("retrying", log)
        self.assertIn("FAILED: step 1b measured no scan-window radius", log)
        self.assertIn(probe_window.EXIT_MEANING[status][:40], log)
        self.assertIn("DependencyNeverSatisfied", log)
        for gone in ("window.txt", "params.resolved.cfg", dp.FALLBACK_RECORD):
            self.assertFalse(os.path.exists(os.path.join(out, gone)), gone)
        w = json.loads(_read(os.path.join(out, "window.json")))
        self.assertEqual(w["exit_status"], status)
        prov = json.loads(_read(os.path.join(out, "search_provenance.json")))
        self.assertNotIn("probe_fallback", prov)
        self.assertEqual(prov["scan_window"]["mode"], "measured", "the plan, untouched")
        return w

    def test_a_handle_that_goes_stale_mid_tail_is_reopened_no_retry_no_traceback(self):
        """The incident itself: the 2nd read of the probe's log raises ESTALE once. Skill
        <= 979dfd5 died here with `OSError: [Errno 116] Stale file handle` at
        `pending += src.read()`; the fix reads its node-local log, reopens on ESTALE at its own
        offset, and measures."""
        with tempfile.TemporaryDirectory() as d:
            out = self._chain(d)
            p = self._run(d, out, "step1b_window.sbatch", estale="transient", FAKE_FIRST="0.5")
            log = p.stdout + p.stderr
            self.assertNotIn("Traceback", log)
            self.assertNotIn("Stale", log)
            self.assertEqual(p.returncode, 0, log)
            self.assertNotIn("retrying", log)
            self.assertEqual(_read(os.path.join(out, "window.txt")).strip(), "7")
            w = json.loads(_read(os.path.join(out, "window.json")))
            self.assertEqual([x["io_error"] for x in w["probes"]], [None] * len(w["probes"]))
            self.assertIsNone(w["failure"])

    def test_estale_on_the_first_attempt_is_retried_and_the_window_measured(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._chain(d)
            p = self._run(d, out, "step1b_window.sbatch", estale="once")
            log = p.stdout + p.stderr
            self.assertEqual(p.returncode, 0, log)
            self.assertNotIn("Traceback", log)
            self.assertIn("retrying once", log)
            self.assertEqual(_read(os.path.join(out, "window.txt")).strip(), "7")
            self.assertFalse(os.path.exists(os.path.join(out, dp.FALLBACK_RECORD)))
            first = json.loads(_read(os.path.join(out, "window.json.attempt1")))
            self.assertEqual(len(first["probes"]), dp.PROBE_MAX_FAILURES)
            self.assertEqual((first["failure"], first["exit_status"]),
                             ("io_error", probe_window.EXIT_IO_ERROR))
            for pr in first["probes"]:
                # the same error reads the same on every run: no per-probe path in it
                self.assertEqual(pr["io_error"], f"OSError: [Errno {errno.ESTALE}] {STALE}")
            # attempt 1's logs are kept beside its evidence, and say what happened
            logs = glob.glob(os.path.join(out, "window_probe.attempt1", "*", "probe.log"))
            self.assertEqual(len(logs), dp.PROBE_MAX_FAILURES, logs)
            self.assertTrue(all("could not be read" in _read(x) for x in logs))
            self.assertEqual(json.loads(_read(os.path.join(out, "window.json")))["window_radius"],
                             7)
            prov = json.loads(_read(os.path.join(out, "search_provenance.json")))
            self.assertEqual(prov["scan_window"]["mode"], "measured")
            self.assertNotIn("probe_fallback", prov)

    def test_estale_every_time_falls_back_and_steps_2_to_5_run_without_a_window(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._chain(d)
            p = self._run(d, out, "step1b_window.sbatch", estale="always")
            log = p.stdout + p.stderr
            self.assertEqual(p.returncode, 0, log)
            self.assertNotIn("Traceback", log)
            self.assertIn("retrying once", log)
            self.assertIn("scan window = auto -- FALLBACK", log)
            self.assertEqual(_read(os.path.join(out, "window.txt")).strip(), "auto")
            fb = json.loads(_read(os.path.join(out, dp.FALLBACK_RECORD)))
            self.assertIn("the probe could not read its own DIA-NN log", fb["reason"])
            self.assertIn("I/O error on the probe log: "
                          f"OSError: [Errno {errno.ESTALE}] {STALE};", fb["reason"])
            prov = json.loads(_read(os.path.join(out, "search_provenance.json")))
            self.assertIn("fallback (probe failed: ", prov["scan_window"]["source"])
            self.assertEqual(prov["scan_window"]["mode"], "fallback_auto")
            self.assertEqual(prov["result"]["scan_window"]["mode"], "fallback_auto")
            resolved = _read(os.path.join(out, "params.resolved.cfg"))
            self.assertNotIn("--window", resolved)
            # steps 2-5: the guard takes `auto` beside the fallback's record, and DIA-NN gets
            # no --window
            for name in ("step2_firstpass.sbatch", "step3_assembly.sbatch",
                         "step4_finalpass.sbatch", "step5_report.sbatch"):
                argv_log = os.path.join(d, f"argv_{name}.txt")
                s = self._run(d, out, name, FAKE_ARGV_LOG=argv_log, FAKE_SLEEP="0",
                              SLURM_ARRAY_TASK_ID="0")
                self.assertNotIn("does not hold", s.stderr, name)
                self.assertNotIn("step 1b", s.stderr, name)
                if os.path.exists(argv_log):
                    self.assertNotIn("--window", _read(argv_log), name)
            self.assertTrue(os.path.exists(os.path.join(d, "argv_step2_firstpass.sbatch.txt")),
                            "step 2 never reached DIA-NN")

            # S5: `auto` without the fallback's record -- deleted, or one that is not a window
            # fallback -- stops steps 2-5 before DIA-NN
            os.rename(os.path.join(out, dp.FALLBACK_RECORD), os.path.join(d, "kept.json"))
            for record in (None, '{"fallback": true, "mass_acc": {"mode": "fallback_default"}}'):
                if record:
                    with open(os.path.join(out, dp.FALLBACK_RECORD), "w") as fh:
                        fh.write(record)
                argv_log = os.path.join(d, "argv_guard.txt")
                s = self._run(d, out, "step2_firstpass.sbatch", FAKE_ARGV_LOG=argv_log,
                              FAKE_SLEEP="0", SLURM_ARRAY_TASK_ID="0")
                self.assertNotEqual(s.returncode, 0)
                self.assertIn("window.txt does not hold a scan-window radius", s.stderr)
                self.assertFalse(os.path.exists(argv_log), "DIA-NN ran on an unrecorded `auto`")

    def test_runs_diann_finished_without_a_radius_fail_step1b_not_fall_back(self):
        """S1: DIA-NN finished the runs and logged no radius -- a wrong FASTA or library, a large
        miscalibration, failed injections. That is the data answering, not the probe failing."""
        with tempfile.TemporaryDirectory() as d:
            out = self._chain(d)
            p = self._run(d, out, "step1b_window.sbatch", FAKE_NO_RADIUS="1")
            w = self._failed_without_fallback(out, p, probe_window.EXIT_NOT_MEASURED)
            self.assertEqual(w["failure"], "not_measured")
            self.assertEqual(len(w["probes"]), dp.PROBE_MAX_FAILURES)

    def test_no_dotnet_fails_step1b_at_once_with_its_own_message(self):
        """S2: DIA-NN cannot read .raw without .NET 8. Every step 2-5 task would fail on the same
        thing, later and noisier -- so step 1b fails, once, and says why."""
        with tempfile.TemporaryDirectory() as d:
            out = self._chain(d)
            p = self._run(d, out, "step1b_window.sbatch", FAKE_NO_DOTNET="1")
            w = self._failed_without_fallback(out, p, probe_window.EXIT_ENVIRONMENT)
            self.assertEqual(w["stopped_because"], "environment")
            self.assertEqual(len(w["probes"]), 1, "no other run can fix the environment")
            self.assertIn("Export DOTNET_ROOT", p.stderr)

    def test_a_refused_mass_accuracy_still_stops_step1b(self):
        """EXIT_REFUSED is never retried or bypassed: a mass accuracy far outside the band means
        the wrong FASTA, species or calibration, and no fallback may hide that."""
        with tempfile.TemporaryDirectory() as d:
            raws = fx._cohort(d)
            for r in raws:
                with open(r + ".log", "w") as fh:
                    fh.write(fx.ms2_as(fx.real_log(os.path.basename(r)[:-4]), 60))
            out, _ = fx.Step1bMassAccChainTests._chain(self, d, raws)
            p = fx.Step1bMassAccChainTests._run_step1b(self, d, out)
            log = p.stdout + p.stderr
            self.assertNotEqual(p.returncode, 0, log)
            self.assertIn("REFUSED", log)
            self.assertNotIn("retrying", log)
            self.assertIn("DependencyNeverSatisfied", log)
            for gone in ("window.txt", "massacc.txt", "params.resolved.cfg",
                         dp.FALLBACK_RECORD):
                self.assertFalse(os.path.exists(os.path.join(out, gone)), gone)


class ProbeFailureTests(unittest.TestCase):
    """probe_window's own classification of a failure, and B1: an implausible mass accuracy
    logged by ANY run is refused -- never fallen back past."""

    @staticmethod
    def _p(name, missing=(), io_error=None, timed_out=False, **vals):
        return dict({"file": f"/x/{name}.raw", "missing": list(missing), "io_error": io_error,
                     "timed_out": timed_out}, **vals)

    def test_the_machinery_falls_back_only_when_every_failed_run_failed_that_way(self):
        io = self._p("a", ["radius"], io_error="OSError: [Errno 116] Stale file handle")
        slow = self._p("b", ["radius"], timed_out=True)
        answered = self._p("c", ["radius"])
        good = self._p("d", radius=7)
        fc = probe_window.failure_class
        self.assertEqual(fc("max_failures", [io, io, io], []), ("io_error", 7))
        self.assertEqual(fc("max_failures", [io, slow, slow], []), ("io_error", 7))
        self.assertEqual(fc("max_failures", [slow, slow, slow], []), ("timed_out", 8))
        self.assertEqual(fc("budget", [good], []), ("timed_out", 8))
        # one run DIA-NN finished without logging it is the data answering
        self.assertEqual(fc("max_failures", [io, io, answered], []), ("not_measured", 6))
        self.assertEqual(fc("no_more_runs", [good], []), ("not_measured", 6))
        self.assertEqual(fc("no_probeable_run", [], []), ("not_measured", 6))
        self.assertEqual(fc("environment", [io], []), ("environment", 4))
        self.assertEqual(fc("max_failures", [io, io, io], ["MS2 60 ppm"]), ("refused", 3))

    def test_implausible_logged_checks_every_run_complete_or_not(self):
        runs = [self._p("a", ["radius"], ms2_ppm=60.0, ms1_ppm=4.0),
                self._p("b", radius=7, ms2_ppm=14.0, ms1_ppm=40.0)]
        bad = probe_window.implausible_logged(runs, ["window", "mass-acc"])
        self.assertEqual(len(bad), 2, bad)
        self.assertIn("MS2 60 ppm (a.raw)", bad[0])
        self.assertIn("MS1 40 ppm (b.raw)", bad[1])
        # a documented level is not measured, so not checked; nor is a window-only probe
        self.assertEqual(len(probe_window.implausible_logged(runs, ["window", "mass-acc"],
                                                              {"ms1_ppm": 7})), 1)
        self.assertEqual(probe_window.implausible_logged(runs, ["window"]), [])

    def _probe(self, d, raws, more=(), env_extra=None, fake=fx.FAKE_DIANN, timeout=60):
        fasta, lib = fx._fasta_lib(d)
        diann = _exe(os.path.join(d, "diann"), fake)
        cfg = os.path.join(d, "resolved.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "probe_window.py"),
                            "--diann", diann, "--raw", *raws, "--fasta", fasta, "--lib", lib,
                            "--threads", "8", "--timeout", str(timeout), "--write-cfg", cfg,
                            "--workdir", os.path.join(d, "w"), *more],
                           capture_output=True, text=True, timeout=240,
                           env=dict(os.environ, DOTNET_ROOT="/opt/fake-dotnet",
                                    **(env_extra or {})))
        return p, cfg

    def test_b1_the_median_run_logs_60_ppm_and_the_others_nothing(self):
        """dda-review B1, scenario 1: the one run that logged everything said 60 ppm; the others
        logged no mass accuracy. Before: too few runs -> exit 1 -> the job fell back and searched
        at 20 ppm with Methods saying "not measured"."""
        with tempfile.TemporaryDirectory() as d:
            raws = fx._cohort(d)
            for r in raws:
                stem = os.path.basename(r)[:-4]
                with open(r + ".log", "w") as fh:
                    fh.write(fx.ms2_as(fx.real_log(stem), 60) if stem == "Ex01162023_10_TT33"
                             else fx.no_ms2(fx.real_log(stem)))
            p, cfg = self._probe(d, raws, ("--measure", "window", "mass-acc", "--ms1-ppm", "7"))
            self.assertEqual(p.returncode, probe_window.EXIT_REFUSED, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(out["failure"], "refused")
            self.assertTrue(any("MS2 60 ppm (Ex01162023_10_TT33.raw)" in x
                                for x in out["mass_acc_refused"]), out["mass_acc_refused"])
            self.assertIn("Refusing to pin", p.stderr)
            self.assertNotIn("--mass-acc", _read(cfg))

    def test_b1_every_run_logs_60_ppm_and_none_a_radius(self):
        """dda-review B1, scenario 2: no run logged everything, so pin_mass_acc never saw the
        60 ppm the runs did log."""
        with tempfile.TemporaryDirectory() as d:
            raws = fx._cohort(d)
            for r in raws:
                log = fx.ms2_as(fx.real_log(os.path.basename(r)[:-4]), 60)
                with open(r + ".log", "w") as fh:
                    fh.write("\n".join(ln for ln in log.splitlines()
                                       if "Scan window radius set to" not in ln) + "\n")
            p, cfg = self._probe(d, raws, ("--measure", "window", "mass-acc", "--ms1-ppm", "7"))
            self.assertEqual(p.returncode, probe_window.EXIT_REFUSED, p.stderr)
            out = json.loads(p.stdout)
            self.assertIsNone(out["window_radius"])
            self.assertEqual([x["ms2_ppm"] for x in out["probes"]], [60.0, 60.0, 60.0])

    def test_a_run_that_goes_silent_is_the_time_limit(self):
        with tempfile.TemporaryDirectory() as d:
            raws = [os.path.join(d, f"r{i}.mzML") for i in range(3)]
            for r in raws:
                with open(r, "wb") as fh:
                    fh.truncate(30 * MB)
            p, _ = self._probe(d, raws, env_extra={"FAKE_HANG": "1"}, fake=FAKE_DIANN, timeout=1)
            self.assertEqual(p.returncode, probe_window.EXIT_TIMED_OUT, p.stderr)
            self.assertEqual(json.loads(p.stdout)["failure"], "timed_out")

    def test_argparse_and_the_probes_own_checks_have_their_own_status(self):
        with tempfile.TemporaryDirectory() as d:
            raws = fx._cohort(d)
            p, _ = self._probe(d, raws, ("--no-such-flag",))
            self.assertEqual(p.returncode, probe_window.EXIT_USAGE, p.stderr)
            p, _ = self._probe(d, raws, ("--ms1-ppm", "7"))       # no --measure mass-acc
            self.assertEqual(p.returncode, probe_window.EXIT_CONFIG, p.stderr)
            self.assertIn("--ms1-ppm/--ms2-ppm pin a documented level", p.stderr)
            p, _ = self._probe(d, raws, ("--", "--dda"))
            self.assertEqual(p.returncode, probe_window.EXIT_CONFIG, p.stderr)

    def test_the_reason_keeps_a_lone_measured_value(self):
        """The fallback discards what one run DID measure; the reason says it."""
        with tempfile.TemporaryDirectory() as d:
            ev = os.path.join(d, "window.json")
            with open(ev, "w") as fh:
                json.dump({"stopped_because": "max_failures", "failure": "io_error",
                           "failed": ["b.raw", "c.raw"],
                           "probes": [self._p("a", radius=7, ms2_ppm=14.0, ms1_ppm=4.1),
                                      self._p("b", ["radius"], io_error="OSError: x"),
                                      self._p("c", ["radius"], io_error="OSError: x")]}, fh)
            why = probe_fallback.reason_from(ev, probe_window.EXIT_IO_ERROR, 2)
            self.assertIn("logged before it failed, not pinned: a.raw radius 7, MS2 14 ppm, "
                          "MS1 4.1 ppm", why)
            self.assertIn("the probe could not read its own DIA-NN log", why)


class ProvenanceModeTests(unittest.TestCase):
    """search_provenance.json `scan_window.mode` / `mass_acc.mode`: STABLE values that FRAN ingests
    as a variable of its DIA-NN vs Spectronaut comparison (fran-5b, 2026-09-30). This class pins
    them: changing a value here is a change to FRAN's data, not a refactor."""

    WINDOW = {"measured", "fallback_auto", "pinned", "auto", "invalid", "unknown"}
    MASS_ACC = {"measured", "fallback_default", "pinned", "pinned_default", "auto", "partial",
                "invalid", "unknown"}

    def test_the_mode_values_are_pinned(self):
        self.assertEqual(set(dp.SCAN_WINDOW_MODES), self.WINDOW)
        self.assertEqual(set(dp.MASS_ACC_MODES), self.MASS_ACC)
        self.assertEqual((probe_fallback.WINDOW_FALLBACK, probe_fallback.MASS_ACC_FALLBACK),
                         ("fallback_auto", "fallback_default"))

    def test_every_mode_is_documented_in_the_provenance_reference(self):
        doc = _read(os.path.join(os.path.dirname(HERE), "references", "environment.md"))
        for key, modes in (("scan_window.mode", self.WINDOW), ("mass_acc.mode", self.MASS_ACC)):
            table = doc.split(f"| `{key}` | meaning |", 1)[1].split("\n\n", 1)[0]
            rows = {ln.split("`")[1] for ln in table.splitlines() if ln.startswith("| `")}
            self.assertEqual(rows, modes, key)

    def _cfg(self, d, body, sidecar=None):
        cfg = os.path.join(d, "p.cfg")
        with open(cfg, "w") as fh:
            fh.write(body)
        if sidecar is not None:
            with open(cfg + ".rationale.json", "w") as fh:
                json.dump(sidecar, fh)
        return cfg

    def test_each_record_says_its_mode(self):
        with tempfile.TemporaryDirectory() as d:
            for body, want in (("--window 7\n", "pinned"), ("--qvalue 0.01\n", "auto"),
                               ("--dda\n", "auto"), ("--window 0\n", "invalid"),
                               ("--window wide\n", "invalid")):
                cfg = self._cfg(d, "--mass-acc 20\n--mass-acc-ms1 7\n" + body)
                self.assertEqual(dp.window_record(dp.mass_acc_status(cfg))["mode"], want, body)
            for body, sidecar, want in (
                    ("--mass-acc 20\n--mass-acc-ms1 7\n", None, "pinned"),
                    ("--qvalue 0.01\n", None, "auto"),
                    ("--mass-acc-ms1 7\n", None, "partial"),
                    ("--mass-acc 0\n--mass-acc-ms1 7\n", None, "invalid"),
                    ("--dda\n--mass-acc 20\n--mass-acc-ms1 10\n",
                     {"mass_accuracy_default": {"--mass-acc": 20}}, "pinned_default")):
                cfg = self._cfg(d, body, sidecar)
                if sidecar is None and os.path.exists(cfg + ".rationale.json"):
                    os.remove(cfg + ".rationale.json")
                rec = dp.mass_acc_record(dp.mass_acc_status(cfg), cfg)
                self.assertEqual(rec["mode"], want, body)
        fb = probe_fallback.fallback(["window", "mass-acc"], {}, "why")
        self.assertEqual((fb["window"]["mode"], fb["mass_acc"]["mode"]),
                         ("fallback_auto", "fallback_default"))

    def test_the_chain_plan_says_measured_and_the_fallback_replaces_it_at_the_top_level(self):
        """A mass-accuracy chain whose probe log goes stale on every run, twice: the documented
        MS1 as given, the SOP MS2 tagged DEFAULT, and steps 2-5 get them with no --window."""
        with tempfile.TemporaryDirectory() as d:
            raws = fx._cohort(d)
            out, info = fx.Step1bMassAccChainTests._chain(self, d, raws)
            self.assertEqual(info["scan_window"]["mode"], "measured")
            self.assertEqual(info["mass_acc"]["mode"], "measured")
            with open(os.path.join(out, "search_provenance.json"), "w") as fh:
                # as run_search.py writes it: scan_window and mass_acc at the top level too
                json.dump({"engine": "diann", "scan_window": info["scan_window"],
                           "mass_acc": info["mass_acc"], "result": info}, fh)
            base = {k: v for k, v in os.environ.items()
                    if k not in ("DOTNET_ROOT", "PROTEOMICS_DOTNET_DIR")}
            p = subprocess.run(["bash", os.path.join(out, "step1b_window.sbatch")], cwd=out,
                               capture_output=True, text=True, timeout=240,
                               env=job_env(d, base=base, **estale_env(d, "always")))
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("mass accuracy = --mass-acc 20 --mass-acc-ms1 7 -- FALLBACK",
                          p.stdout + p.stderr)
            prov = json.loads(_read(os.path.join(out, "search_provenance.json")))
            self.assertEqual(prov["scan_window"]["mode"], "fallback_auto")
            self.assertEqual(prov["mass_acc"]["mode"], "fallback_default")
            self.assertEqual(prov["result"]["mass_acc"], prov["mass_acc"])
            self.assertEqual(prov["mass_acc"]["default"], {"--mass-acc": 20.0})
            self.assertEqual(prov["mass_acc"]["documented"], {"--mass-acc-ms1": 7.0})
            self.assertEqual(_read(os.path.join(out, "massacc.txt")).strip(),
                             "--mass-acc 20 --mass-acc-ms1 7")
            p2, argv_log = fx.Step1bMassAccChainTests._run_step(self, d, out,
                                                                "step2_firstpass.sbatch")
            self.assertNotIn("step 1b", p2.stderr)
            call = _read(argv_log)
            self.assertIn("--mass-acc 20 --mass-acc-ms1 7", call)
            self.assertNotIn("--window", call)


class CautionTests(unittest.TestCase):
    """A search that fell back is a CAUTION wherever the results are read, not only in the
    Methods: AUDIT.md (and so the report's Audit & caveats), the report's own Results-at-a-glance
    callout, and the run record's Data Quality Notes -- one wording, probe_fallback.caution()."""

    REASON = ("the runs tried did not log it; I/O error on the probe log: OSError: [Errno 116] "
              "Stale file handle; 2 attempts, last exit 1")

    def _prov(self, d, window="fallback_auto", mass_acc="fallback_default"):
        path = os.path.join(d, "search", "search_provenance.json")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            json.dump({"engine": "diann", "search_mode": "single_shot",
                       "scan_window": {"mode": window},
                       "mass_acc": {"mode": mass_acc, "pin_as": "--mass-acc 20 --mass-acc-ms1 7"},
                       **({"probe_fallback": {"fallback": True, "reason": self.REASON}}
                          if "fallback" in window + mass_acc else {})}, fh)
        return path

    def test_the_caution_wording(self):
        text = probe_fallback.caution("fallback_auto", "measured", "x")
        self.assertTrue(text.startswith("CAUTION: the DIA-NN scan window was NOT measured"))
        self.assertIn("chose the window itself for each run", text)
        self.assertIn("advises against combining", text)
        self.assertIn("Re-run the search if the cause was transient", text)
        self.assertNotIn("mass accuracy", text)
        both = probe_fallback.caution("fallback_auto", "fallback_default", "x", "--mass-acc 20")
        self.assertIn("the mass accuracy was NOT measured", both)
        self.assertIn("(--mass-acc 20)", both)
        for sw, ma in (("measured", "measured"), ("auto", "pinned"), (None, None)):
            self.assertIsNone(probe_fallback.caution(sw, ma, "x"))

    def test_one_reader_for_the_methods_and_the_caution_old_provenance_included(self):
        """M2 (rule 3): the Methods sentence and every CAUTION key off the stable mode, through
        one reader (probe_fallback.fallback_modes). A provenance from 7f907ed / 8290ab0 has
        `probe_fallback` and no mode (the first HIVE test chains): the reader derives
        fallback_auto / fallback_default from it, so the Methods and the CAUTION agree."""
        old = {"engine": "diann", "scan_window": {"source": "fallback (probe failed: x)"},
               "probe_fallback": {"fallback": True, "reason": self.REASON,
                                  "window": {"value": "auto", "fallback": True},
                                  "mass_acc": {"pin_as": "--mass-acc 20 --mass-acc-ms1 7"}},
               "result": {"mass_acc": {"fixed": True, "measured": False}}}
        self.assertEqual(probe_fallback.fallback_modes(old), ("fallback_auto", "fallback_default"))
        new = {"scan_window": {"mode": "measured"}, "mass_acc": {"mode": "pinned"}}
        self.assertEqual(probe_fallback.fallback_modes(new), ("measured", "pinned"))
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "search", "search_provenance.json")
            os.makedirs(os.path.dirname(path))
            with open(path, "w") as fh:
                json.dump(old, fh)
            fb = make_methods.probe_fallback_record(old)
            self.assertEqual((fb["window"], fb["mass_acc"]), (True, True))
            self.assertIn("scan window", make_methods.probe_fallback_sentence({"probe_fallback": fb}))
            md, js = self._audit(d, path)
            self.assertIn("CAUTION: the DIA-NN scan window was NOT measured", md)
            self.assertIn("(--mass-acc 20 --mass-acc-ms1 7)", md)
            note = make_analysis_html.search_measurement_note(
                {"input": os.path.join(os.path.dirname(path), "report.parquet")})
            self.assertIn("CAUTION: the DIA-NN scan window", note["text"])
            params = {"scan_window": old["scan_window"], "mass_acc": None,
                      "probe_fallback": old["probe_fallback"]}
            self.assertIn("scan window was NOT measured", record_run.measurement_caution(params))
        # a provenance with neither: nobody says anything
        self.assertIsNone(make_methods.probe_fallback_record({"engine": "diann"}))
        self.assertIsNone(probe_fallback.provenance_caution({"engine": "diann"}))

    def test_a_malformed_record_never_stops_the_report_or_the_audit(self):
        """M3: scan_window / mass_acc / probe_fallback that are not objects."""
        for bad in ({"scan_window": "x", "mass_acc": 7, "probe_fallback": "y"},
                    {"scan_window": ["fallback_auto"], "result": "z"},
                    {"scan_window": {"mode": 3}, "probe_fallback": {"window": "auto"}}):
            with self.subTest(bad=bad), tempfile.TemporaryDirectory() as d:
                self.assertEqual(probe_fallback.fallback_modes(bad), (None, None))
                self.assertIsNone(probe_fallback.provenance_caution(bad))
                path = os.path.join(d, "search", "search_provenance.json")
                os.makedirs(os.path.dirname(path))
                with open(path, "w") as fh:
                    json.dump(bad, fh)
                self.assertIsNone(make_analysis_html.search_measurement_note(
                    {"input": os.path.join(os.path.dirname(path), "report.parquet")}))
                md, js = self._audit(d, path)
                self.assertFalse([x for x in js["findings"]
                                  if x["check"] == "search_measurement"])
                self.assertIsNone(make_methods.probe_fallback_record(bad))

    def _audit(self, d, prov):
        out = os.path.join(d, "AUDIT.md")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"),
                            "--out", out, "--search-prov", prov],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        return _read(out), json.loads(_read(os.path.join(d, "AUDIT.json")))

    def test_audit_md_carries_the_caution_as_a_warning(self):
        with tempfile.TemporaryDirectory() as d:
            md, js = self._audit(d, self._prov(d))
            f = [x for x in js["findings"] if x["check"] == "search_measurement"]
            self.assertEqual(len(f), 1, js["findings"])
            self.assertEqual(f[0]["status"], "WARN")
            self.assertEqual(f[0]["detail"]["scan_window_mode"], "fallback_auto")
            self.assertIn("**search_measurement** — CAUTION: the DIA-NN scan window was NOT "
                          "measured", md)
            self.assertIn("Stale file handle", md)
        with tempfile.TemporaryDirectory() as d:
            md, js = self._audit(d, self._prov(d, "measured", "measured"))
            self.assertFalse([x for x in js["findings"] if x["check"] == "search_measurement"])

    def test_the_report_shows_it_at_a_glance_whoever_wrote_the_report(self):
        with tempfile.TemporaryDirectory() as d:
            prov = self._prov(d)
            note = make_analysis_html.search_measurement_note(
                {"input": os.path.join(os.path.dirname(prov), "report.parquet")})
            self.assertEqual(note["kind"], "warning")
            self.assertTrue(note["text"].startswith("CAUTION: the DIA-NN scan window"))
            # the page itself, with no AI report and no DE tables
            tables = os.path.join(d, "tables")
            os.makedirs(tables)
            with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
                json.dump({"input": os.path.join(os.path.dirname(prov), "report.parquet")}, fh)
            html = os.path.join(d, "Analysis_Report.html")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                                "--tables", tables, "--out", html, "--no-pdf"],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(p.returncode, 0, p.stderr)
            for f in (html, os.path.join(d, "Analysis_Report.md")):
                self.assertIn("CAUTION: the DIA-NN scan window was NOT measured", _read(f), f)
        with tempfile.TemporaryDirectory() as d:
            prov = self._prov(d, "measured", "measured")
            self.assertIsNone(make_analysis_html.search_measurement_note(
                {"input": os.path.join(os.path.dirname(prov), "report.parquet")}))

    def _fallen_back(self, d):
        """A search output directory whose probe fell back, with its record."""
        out = os.path.join(d, "search")
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, dp.FALLBACK_RECORD), "w") as fh:
            json.dump(probe_fallback.fallback(["window", "mass-acc"], {"ms1_ppm": 7.0},
                                              self.REASON), fh, indent=2)
        return out

    def test_the_status_views_say_it(self):
        """checkpoint.py status, the Slack post and watch_run.sh --all -- where a resumed
        session and the Core channel see a search's state."""
        import notify_slack
        with tempfile.TemporaryDirectory() as d:
            out = self._fallen_back(d)
            self.assertTrue(probe_fallback.caution_for(out).startswith("CAUTION: the DIA-NN "
                                                                       "scan window"))
            env = job_env(d)
            sess = os.path.join(d, "session")
            os.makedirs(sess)
            ck = os.path.join(SCRIPTS, "checkpoint.py")
            r = subprocess.run([sys.executable, ck, "record", "--session", sess, "--stage",
                                "search", "--jobs", "1", "--report",
                                os.path.join(out, "report.parquet")],
                               capture_output=True, text=True, env=env, timeout=60)
            self.assertEqual(r.returncode, 0, r.stderr)
            r = subprocess.run([sys.executable, ck, "status", "--session", sess],
                               capture_output=True, text=True, env=env, timeout=60)
            st = json.loads(r.stdout)
            self.assertTrue(any(w.startswith("CAUTION: the DIA-NN scan window was NOT measured")
                                for w in st.get("warnings") or []), st)
            # Slack: the headline and a block of its own
            f = notify_slack.search_facts(out, exit_code=0)
            self.assertTrue(f["caution"].startswith("CAUTION: the DIA-NN scan window"))
            msg = notify_slack.render(f)
            self.assertIn("CAUTION: fell back, not measured", msg["text"])
            self.assertTrue(any("scan window was NOT measured" in json.dumps(b)
                                for b in msg["blocks"]))
            with open(os.path.join(out, "jobs.txt"), "w") as fh:
                fh.write("1\n")
            w = subprocess.run(["bash", os.path.join(SCRIPTS, "watch_run.sh"), "--all", out],
                               capture_output=True, text=True, env=env, timeout=60)
            js = json.loads(w.stdout)
            self.assertTrue(js["probe_fallback"])
            self.assertIn("FELL BACK", js["caution"])
        with tempfile.TemporaryDirectory() as d:          # nothing fell back: nothing said
            out = os.path.join(d, "search")
            os.makedirs(out)
            self.assertIsNone(probe_fallback.caution_for(out))
            self.assertIsNone(notify_slack.search_facts(out, exit_code=0)["caution"])
            with open(os.path.join(out, "jobs.txt"), "w") as fh:
                fh.write("1\n")
            w = subprocess.run(["bash", os.path.join(SCRIPTS, "watch_run.sh"), "--all", out],
                               capture_output=True, text=True, env=job_env(d), timeout=60)
            self.assertNotIn("probe_fallback", json.loads(w.stdout))

    def test_a_stale_fallback_record_never_describes_the_next_search(self):
        """dda-review N1: a pinned-window chain generated into a folder holding an earlier
        chain's probe_fallback.json. Generation removes every probe output of the earlier
        search, and the status views read the provenance before any record file."""
        with tempfile.TemporaryDirectory() as d:
            runs = []
            for i in range(3):
                p = os.path.join(d, f"run{i}.mzML")
                with open(p, "wb") as fh:
                    fh.truncate((30 + i) * MB)
                runs.append(p)
            fasta, lib = fx._fasta_lib(d)
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n--mass-acc 20\n--mass-acc-ms1 7\n--window 7\n")
            out = self._fallen_back(d)                      # the earlier search's record
            for name in ("window.txt", "massacc.txt", "window.json", "window.json.attempt1"):
                with open(os.path.join(out, name), "w") as fh:
                    fh.write("from an earlier search\n")
            self.assertTrue(probe_fallback.caution_for(out), "fixture")
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", _exe(os.path.join(d, "diann"), FAKE_DIANN),
                                "--raw", *runs, "--fasta", fasta, "--out", out, "--cfg", cfg,
                                "--threads-per-file", "4"],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(g.returncode, 0, g.stderr)
            self.assertIn("do not describe this search: set aside as", g.stderr)
            for name in ("window.txt", "massacc.txt", "window.json", "window.json.attempt1",
                         dp.FALLBACK_RECORD):
                self.assertFalse(os.path.exists(os.path.join(out, name)), name)
                # set aside, never deleted: that search's evidence
                self.assertEqual(len(glob.glob(os.path.join(out, name + ".stale-*"))), 1, name)
            info = json.loads(g.stdout)
            self.assertEqual(info["scan_window"]["mode"], "pinned")
            self.assertIsNone(probe_fallback.caution_for(out))
            # and with a stale record back beside this search's provenance, the provenance wins
            with open(os.path.join(out, "search_provenance.json"), "w") as fh:
                json.dump({"engine": "diann", "scan_window": info["scan_window"],
                           "mass_acc": info["mass_acc"], "result": info}, fh, indent=2)
            self._fallen_back(d)
            self.assertIsNone(probe_fallback.caution_for(out))
            with open(os.path.join(out, "jobs.txt"), "w") as fh:
                fh.write("1\n")
            w = subprocess.run(["bash", os.path.join(SCRIPTS, "watch_run.sh"), "--all", out],
                               capture_output=True, text=True, env=job_env(d), timeout=60)
            self.assertNotIn("probe_fallback", json.loads(w.stdout))

    def test_the_run_record_copies_the_fallback_record_only_when_the_provenance_says_so(self):
        """dda-review N1: record_run.py copied whatever probe_fallback.json sat in the folder into
        the run record. It now asks the provenance (fallback_modes, the one reader) -- a new
        provenance's mode, or an old one's probe_fallback -- and copies the file only then."""
        import types
        a = types.SimpleNamespace(issues_tag=None, skill_version="test", file_cap_mb=50,
                                  zip_cap_gb=1, detect_json=None, exit_code=0, json=None,
                                  meta=None, status="completed")
        cases = (({"scan_window": {"mode": "pinned"}, "mass_acc": {"mode": "pinned"}}, False),
                 ({"scan_window": {"mode": "fallback_auto"}}, True),
                 ({"mass_acc": {"mode": "fallback_default"}}, True),
                 ({"scan_window": {"source": "x"},                     # before the modes
                   "probe_fallback": {"window": {"value": "auto"}}}, True),
                 (None, False))                                         # no provenance at all
        for prov, copied in cases:
            with self.subTest(prov=prov), tempfile.TemporaryDirectory() as d:
                with open(os.path.join(d, dp.FALLBACK_RECORD), "w") as fh:
                    fh.write("{}\n")
                if prov is not None:
                    with open(os.path.join(d, "search_provenance.json"), "w") as fh:
                        json.dump(dict(prov, engine="diann"), fh, indent=2)
                plan = record_run.new_plan("search-done", a)
                with mock.patch.dict(os.environ, {"RECORD_RUN": "off", "SKILL_RUNS_DIR": "off"}):
                    record_run.plan_search(plan, d, a, time.time() + 60)
                rels = {c["rel"] for c in plan["copies"]}
                self.assertEqual("output/search/probe_fallback.json" in rels, copied, rels)
        # watch_run.sh --all reads the provenance the same way, old ones included
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "jobs.txt"), "w") as fh:
                fh.write("1\n")
            with open(os.path.join(d, "search_provenance.json"), "w") as fh:
                json.dump(cases[3][0], fh, indent=2)
            w = subprocess.run(["bash", os.path.join(SCRIPTS, "watch_run.sh"), "--all", d],
                               capture_output=True, text=True, env=job_env(d), timeout=60)
            self.assertTrue(json.loads(w.stdout).get("probe_fallback"), w.stdout)

    def test_a_refused_generation_leaves_a_completed_search_alone(self):
        """dda-review R1: the earlier search's probe outputs were cleared BEFORE the generator's
        refusals, so a REFUSED generation deleted a completed search's window.txt, massacc.txt
        and window.json and left its provenance pointing at files that were gone."""
        with tempfile.TemporaryDirectory() as d:
            runs = []
            for i in range(3):
                p = os.path.join(d, f"run{i}.mzML")
                with open(p, "wb") as fh:
                    fh.truncate((30 + i) * MB)
                runs.append(p)
            fasta, _ = fx._fasta_lib(d)
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n")                    # mass accuracy not pinned: refused
            out = os.path.join(d, "out")
            os.makedirs(out)
            kept = {"window.txt": "7\n", "massacc.txt": "--mass-acc 20 --mass-acc-ms1 7\n",
                    "window.json": '{"window_radius": 7}\n', "report.parquet": "PAR1"}
            for name, text in kept.items():
                with open(os.path.join(out, name), "w") as fh:
                    fh.write(text)
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", _exe(os.path.join(d, "diann"), FAKE_DIANN),
                                "--raw", *runs, "--fasta", fasta, "--out", out, "--cfg", cfg,
                                "--threads-per-file", "4"],
                               capture_output=True, text=True, timeout=120)
            self.assertNotEqual(g.returncode, 0, "fixture: this cfg is refused")
            self.assertNotIn("set aside", g.stderr)
            for name, text in kept.items():
                self.assertEqual(_read(os.path.join(out, name)), text, name)
            self.assertEqual(glob.glob(os.path.join(out, "*.stale-*")), [])

    def test_the_run_record_says_it_once(self):
        params = {"scan_window": {"mode": "fallback_auto"}, "mass_acc": {"mode": "measured"},
                  "probe_fallback": {"reason": self.REASON}}
        rec = {"search": {"out_dir": "/x/out", "parameters": params, "status": "completed"},
               "prot": {"prot": "PROT_0001"}}
        notes = record_run.data_quality_notes(rec)
        hits = [n for n in notes if "scan window was NOT measured" in n["what"]]
        self.assertEqual(len(hits), 1, notes)
        self.assertEqual(hits[0]["severity"], "WARNING")
        self.assertIn("Stale file handle", hits[0]["what"])
        self.assertIn("## Data Quality Notes", "\n".join(record_run.render_dq(notes)))
        # AUDIT.md already said it: the audit's finding is the one note
        rec["analysis"] = {"audit": [{"check": "search_measurement", "status": "WARN",
                                      "message": probe_fallback.caution(
                                          "fallback_auto", "measured", "x")}]}
        notes = record_run.data_quality_notes(rec)
        self.assertEqual(len([n for n in notes if "scan window was NOT measured" in n["what"]]), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
