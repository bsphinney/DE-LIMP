#!/usr/bin/env python3
"""
The single-shot DIA-NN search (<= 5 files, or no SLURM chain) must be as safe to submit as
the 5-step chain. Found on the FRAN DIA-NN 2.7.0 pilot on HIVE (2026-09-16), where every
one of these had to be caught by hand before submission:

  1. QUEUE. `--partition low --account publicgrp --qos publicgrp-low-qos` was passed to
     run_search.py, but run_diann() called emit_sbatch() without them, so slurm_queue()
     re-detected and wrote genome-center-grp/high into both job headers. The headers were
     hand-edited. --requeue was keyed on `partition == "low"` only, while the chain's
     header() also requeues on a public QOS -- two rules for one concept.
  2. OUTPUT. DIA-NN exits 0 on a fatal error (references/diann_parallel.md). The chain
     asserts every artefact; the two-job single-shot path asserted nothing, so a library
     job that wrote no library released the search job, and a run DIA-NN could not read
     dropped out of report.parquet while SLURM said COMPLETED.
  3. --rt-profiling. It was stripped from the search job, which with --reanalyse is the job
     that builds the empirical library the flag configures.

How the cfg is spliced into the search job (quoting, comments, globs) is main's
diann_parallel.cfg_tokens/bash_flags, tested in tests/test_cfg_reader_quoting.py; it is not
re-tested here.

Stdlib only, like the rest of the suite: CI installs no pyarrow. The report check has a
pyarrow-free route (DIA-NN's own .stats.tsv) precisely so it can run -- and be tested --
without it; the pyarrow route is tested only where pyarrow exists.
"""
import argparse
import io
import json
import os
import re
import shlex
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import diann_parallel  # noqa: E402
import run_search  # noqa: E402

try:
    import pyarrow  # noqa: F401
    HAVE_PYARROW = True
except ImportError:
    HAVE_PYARROW = False

LIBFREE_CFG = [
    "--qvalue 0.01", "--matrices", "--xic 10", "--mobilograms",
    "--fasta-search", "--gen-spec-lib", "--predictor", "--reanalyse", "--rt-profiling",
    "--cut K*,R*", "--var-mod UniMod:1,42.010565,*n",
    "--mass-acc 23.3", "--mass-acc-ms1 10",
]

# A stand-in for DIA-NN that behaves the way DIA-NN does on a fatal error: it exits 0 no
# matter what. It writes only what the FAKE_* environment asks for, and records its argv.
# The stats file mirrors DIA-NN 2.7.0's <report>.stats.tsv: one row per run, File.Name is
# the input path as given, Precursors.Identified is 0 for a run that produced nothing.
FAKE_DIANN = r'''#!/bin/bash
[ -n "${FAKE_ARGV:-}" ] && printf '%s\n' "$@" > "$FAKE_ARGV"
out_lib=""; out=""; files=()
while [ $# -gt 0 ]; do
  case "$1" in
    --out-lib) out_lib="$2"; shift ;;
    --out) out="$2"; shift ;;
    --f) files+=("$2"); shift ;;
  esac
  shift
done
if [ -n "$out_lib" ] && [ "${FAKE_WRITE_LIB:-0}" = 1 ]; then
  printf 'predicted library bytes' > "${out_lib}.predicted.speclib"
fi
if [ -n "$out" ] && [ "${FAKE_EMPTY_REPORT:-0}" = 1 ]; then
  : > "$out"                      # created the report, then died: 0 bytes, still exit 0
  exit 0
fi
if [ -n "$out" ] && [ -n "${FAKE_REPORT_RUNS+x}" ]; then
  "$FAKE_PY" - "$out" "$FAKE_REPORT_RUNS" "${files[@]}" <<'PY'
import os, sys
out, keep, inputs = sys.argv[1], [r for r in sys.argv[2].split(",") if r], sys.argv[3:]
try:
    import pyarrow as pa, pyarrow.parquet as pq
    pq.write_table(pa.table({"Run": keep, "Protein.Group": ["P1"] * len(keep)}), out)
except ImportError:
    open(out, "wb").write(b"PAR1 placeholder: no pyarrow in this interpreter")
stem = out[:-len(".parquet")] if out.endswith(".parquet") else out
if os.environ.get("FAKE_NO_STATS") == "1":      # --no-stats in the cfg
    sys.exit(0)
with open(stem + ".stats.tsv", "w") as fh:
    fh.write("File.Name\tPrecursors.Identified\tProteins.Identified\n")
    for f in inputs:
        name = os.path.splitext(os.path.basename(f.rstrip("/")))[0]
        n = 1500 if name in keep else 0
        fh.write(f"{f}\t{n}\t{n // 10}\n")
PY
fi
exit 0
'''


def _read(path):
    with open(path) as fh:
        return fh.read()


def _write(path, text, mode=None):
    with open(path, "w") as fh:
        fh.write(text)
    if mode:
        os.chmod(path, mode)
    return path


class _Workspace:
    """A throwaway search: fake DIA-NN, tools.json, bundle, cfg, FASTA and input runs."""

    def __init__(self, d, cfg_lines=LIBFREE_CFG, n_runs=2, ext=".mzML"):
        self.d = d
        self.diann = _write(os.path.join(d, "fake-diann"), FAKE_DIANN, 0o755)
        self.tools = _write(os.path.join(d, "tools.json"),
                            json.dumps({"diann": self.diann, "versions": {"diann": "2.7.0"}}))
        self.bundle = _write(os.path.join(d, "manifest.json"),
                             json.dumps({"engine": {"name": "diann"}, "acquisition": "DIA"}))
        self.cfg = _write(os.path.join(d, "diann.cfg"), "\n".join(cfg_lines) + "\n")
        self.fasta = _write(os.path.join(d, "search.fasta"), ">sp|P1|X\nPEPTIDER\n")
        # mzML, not .raw: .raw makes run_search provision .NET, which needs the network.
        self.runs = [_write(os.path.join(d, f"sample_{i}{ext}"), "") for i in range(n_runs)]
        self.out = os.path.join(d, "search_out")
        self.job = os.path.join(d, "diann_job.sh")

    def generate(self, *queue, check=True):
        p = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "run_search.py"),
             "--tools", self.tools, "--bundle", self.bundle, "--params", self.cfg,
             "--fasta", self.fasta, "--out", self.out, "--files", *self.runs,
             "--threads", "4", "--sbatch", self.job, *queue],
            capture_output=True, text=True)
        if check and p.returncode != 0:
            raise AssertionError(f"run_search.py failed:\n{p.stdout}\n{p.stderr}")
        return p

    @property
    def lib_job(self):
        return os.path.join(self.d, "diann_job_1_lib.sh")

    @property
    def search_job(self):
        return os.path.join(self.d, "diann_job_2_search.sh")

    def run_job(self, script, **env):
        e = dict(os.environ, FAKE_PY=sys.executable)
        e.update({k: str(v) for k, v in env.items()})
        return subprocess.run(["bash", script], capture_output=True, text=True, env=e)

    def run_inline(self, **env):
        """run_search.py with no --sbatch: the search runs in run_search.py's own shell."""
        e = dict(os.environ, FAKE_PY=sys.executable)
        e.update({k: str(v) for k, v in env.items()})
        return subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "run_search.py"),
             "--tools", self.tools, "--bundle", self.bundle, "--params", self.cfg,
             "--fasta", self.fasta, "--out", self.out, "--files", *self.runs,
             "--threads", "4", "--allow-inline"],
            capture_output=True, text=True, env=e)


def _headers(path):
    return [l for l in _read(path).splitlines() if l.startswith("#SBATCH")]


LOW = ("--partition", "low", "--account", "publicgrp", "--qos", "publicgrp-low-qos")
HIGH = ("--partition", "high", "--account", "genome-center-grp")


class QueueTests(unittest.TestCase):
    def test_the_two_job_chain_carries_the_queue_the_user_passed(self):
        """The pilot's exact invocation. Each header line exactly once, in both jobs."""
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*LOW)
            for job in (w.lib_job, w.search_job):
                h = _headers(job)
                for want in ("#SBATCH --partition=low", "#SBATCH --account=publicgrp",
                             "#SBATCH --qos=publicgrp-low-qos", "#SBATCH --requeue"):
                    self.assertEqual(h.count(want), 1, f"{os.path.basename(job)}: {want}\n{h}")
                other = [l for l in h if re.match(r"#SBATCH --(partition|account|qos)=", l)
                         and l.split("=", 1)[1] not in ("low", "publicgrp", "publicgrp-low-qos")]
                self.assertEqual(other, [], os.path.basename(job))

    def test_the_priority_queue_is_not_marked_preemptible(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*HIGH)
            for job in (w.lib_job, w.search_job):
                h = _headers(job)
                self.assertIn("#SBATCH --partition=high", h)
                self.assertIn("#SBATCH --account=genome-center-grp", h)
                self.assertNotIn("#SBATCH --requeue", h)
                # high needs no --qos (measured, see test_hive_submission_guards); none invented
                self.assertFalse([l for l in h if l.startswith("#SBATCH --qos")], h)

    def test_the_single_job_path_carries_the_queue_too(self):
        """A cfg that is not library-free (a --lib search, or --one-step) takes the one-job
        branch, which had the same missing arguments."""
        cfg = [l for l in LIBFREE_CFG if l not in ("--fasta-search", "--predictor")]
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=cfg)
            w.generate(*LOW)
            self.assertFalse(os.path.exists(w.lib_job))
            h = _headers(w.job)
            self.assertIn("#SBATCH --partition=low", h)
            self.assertIn("#SBATCH --account=publicgrp", h)
            self.assertIn("#SBATCH --qos=publicgrp-low-qos", h)
            self.assertIn("#SBATCH --requeue", h)

    def test_requeue_follows_the_same_rule_as_the_chain(self):
        """One definition: a single-shot job and a chain step on the same queue must agree
        on whether a preemption requeues it or loses it."""
        cases = [("low", None), ("low", "publicgrp-low-qos"), ("high", None),
                 ("high", "genome-center-grp-high-qos"), ("high", "publicgrp-high-qos")]
        with tempfile.TemporaryDirectory() as d:
            for part, qos in cases:
                path = os.path.join(d, f"{part}_{qos}.sh")
                run_search.emit_sbatch(path, "true", d, 1, job="x", partition=part,
                                       account="acct", qos=qos)
                single = "#SBATCH --requeue" in _headers(path)
                chain = "#SBATCH --requeue" in diann_parallel.header(
                    "x", 1, 1, 1, part, "acct", qos=qos).splitlines()
                self.assertEqual(single, chain, (part, qos))

    def test_every_emit_sbatch_call_forwards_the_queue(self):
        """Sage, AlphaDIA, FragPipe and single-file Radiant emit through the same function and
        had the same omission; a new engine must not reintroduce it."""
        src = _read(os.path.join(SCRIPTS, "run_search.py"))
        calls = []
        for m in re.finditer(r"^[ \t]*emit_sbatch\(", src, re.M):   # statements, not prose
            depth, i = 0, m.end() - 1
            for i in range(m.end() - 1, len(src)):          # balanced parens: args may nest
                depth += {"(": 1, ")": -1}.get(src[i], 0)
                if depth == 0:
                    break
            calls.append(src[m.start():i + 1])
        self.assertGreaterEqual(len(calls), 5, calls)
        for c in calls:
            self.assertRegex(c, r"\*\*\(?queue\b", f"emit_sbatch call without the queue: {c}")


# HIVE's real associations for a facility member who is also in publicgrp, read with
# `sacctmgr -nP show assoc user=brettsp format=account,partition,qos` on 2026-09-16.
HIVE_ASSOCIATIONS = """genome-center-grp|gpu-a100|genome-center-grp-gpu-a100-qos
genome-center-grp|high|genome-center-grp-high-qos
publicgrp|high|publicgrp-high-qos
publicgrp|low|publicgrp-low-qos
"""


class PartialQueueTests(unittest.TestCase):
    """Forwarding --partition/--account/--qos to every job made a PARTIAL override dangerous.

    slurm_queue() filled the missing fields from the first preferred association, not one
    that matched what was given. The review generated these against HIVE's associations and
    checked each header with `srun --test-only`:
      --partition low  -> low / genome-center-grp / genome-center-grp-high-qos
                          "Invalid account or account/partition combination specified"
      --qos publicgrp-low-qos -> high / genome-center-grp / publicgrp-low-qos
                          "Invalid qos specification"
    Before the flags were forwarded these paths ignored them and wrote a valid detected
    queue; the 5-step chain had the same fill-in and still does unless this is fixed there."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        fakebin = os.path.join(self._tmp.name, "bin")
        os.makedirs(fakebin)
        _write(os.path.join(fakebin, "sacctmgr"),
               "#!/bin/sh\ncat <<'ASSOC'\n" + HIVE_ASSOCIATIONS + "ASSOC\n", 0o755)
        env = mock.patch.dict(os.environ, {"PATH": fakebin + os.pathsep + os.environ["PATH"],
                                           "USER": "tester"})
        env.start()
        self.addCleanup(env.stop)
        self.addCleanup(self._tmp.cleanup)

    def _queue(self, path):
        h = _headers(path)
        get = lambda k: next((l.split("=", 1)[1] for l in h if l.startswith(f"#SBATCH --{k}=")),
                             None)
        return get("partition"), get("account"), get("qos"), "#SBATCH --requeue" in h

    def _both_jobs(self, *flags):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*flags)
            got = {self._queue(w.lib_job), self._queue(w.search_job)}
            self.assertEqual(len(got), 1, got)
            return got.pop()

    def test_a_partition_alone_takes_the_account_that_has_it(self):
        self.assertEqual(self._both_jobs("--partition", "low"),
                         ("low", "publicgrp", "publicgrp-low-qos", True))

    def test_a_qos_alone_takes_the_association_it_belongs_to(self):
        self.assertEqual(self._both_jobs("--qos", "publicgrp-low-qos"),
                         ("low", "publicgrp", "publicgrp-low-qos", True))

    def test_an_account_alone_stays_on_that_account(self):
        part, acct, qos = run_search.slurm_queue(account="publicgrp")
        self.assertEqual(acct, "publicgrp")
        self.assertIn(f"{acct}|{part}|{qos}", HIVE_ASSOCIATIONS)
        part, acct, qos = run_search.slurm_queue(account="genome-center-grp")
        self.assertEqual((part, acct), ("high", "genome-center-grp"))

    def test_a_combination_no_association_has_fails_loudly(self):
        """Nothing can be filled in to make this valid, so no job is written at all."""
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            p = w.generate("--partition", "low", "--qos", "genome-center-grp-high-qos",
                           check=False)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("publicgrp|low", p.stderr)          # names what the user CAN use
            self.assertFalse(os.path.exists(w.lib_job))
            self.assertFalse(os.path.exists(w.search_job))

    def test_partition_and_account_get_the_qos_the_chain_gives_them(self):
        """diann_parallel adds publicgrp-low-qos for low/publicgrp; emit_sbatch did not, so the
        same flags wrote different headers on the two paths."""
        self.assertEqual(self._both_jobs("--partition", "low", "--account", "publicgrp"),
                         ("low", "publicgrp", "publicgrp-low-qos", True))

    def test_the_chain_resolves_a_partial_queue_the_same_way(self):
        for given in ({"partition": "low"}, {"qos": "publicgrp-low-qos"},
                      {"partition": "low", "account": "publicgrp"}):
            with self.subTest(**given), tempfile.TemporaryDirectory() as d:
                w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], n_runs=3, ext=".d")
                a = argparse.Namespace(partition=None, account=None, qos=None,
                                       max_simultaneous=None)
                for k, v in given.items():
                    setattr(a, k, v)
                os.makedirs(w.out, exist_ok=True)
                run_search.run_diann_parallel(w.diann, w.cfg, w.runs, w.fasta, w.out, 16, a)
                for step in ("step1_libpred.sbatch", "step2_firstpass.sbatch",
                             "step5_report.sbatch"):
                    self.assertEqual(self._queue(os.path.join(w.out, step)),
                                     ("low", "publicgrp", "publicgrp-low-qos", True), step)


class LibraryJobGuardTests(unittest.TestCase):
    def test_a_library_job_that_writes_no_library_fails(self):
        """Otherwise afterok releases a search against a library that does not exist."""
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*HIGH)
            p = w.run_job(w.lib_job, FAKE_WRITE_LIB=0)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("diann_lib.predicted.speclib", p.stderr)

    def test_a_library_job_that_writes_the_library_passes(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*HIGH)
            p = w.run_job(w.lib_job, FAKE_WRITE_LIB=1)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertTrue(os.path.exists(os.path.join(w.out, "diann_lib.predicted.speclib")))


class SearchJobGuardTests(unittest.TestCase):
    def test_no_report_fails(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*HIGH)
            p = w.run_job(w.search_job)                    # FAKE_REPORT_RUNS unset: no report
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("report.parquet", p.stderr)

    def test_a_run_missing_from_the_report_fails_and_is_named(self):
        """The silent failure: DIA-NN could not read sample_1, left it out, exited 0."""
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*HIGH)
            p = w.run_job(w.search_job, FAKE_REPORT_RUNS="sample_0")
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("sample_1", p.stdout + p.stderr)
            self.assertIn("1 of 2", p.stdout + p.stderr)

    def test_every_run_present_passes(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*HIGH)
            p = w.run_job(w.search_job, FAKE_REPORT_RUNS="sample_0,sample_1")
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("all 2", p.stdout)

    def test_the_single_job_path_checks_the_report_too(self):
        cfg = [l for l in LIBFREE_CFG if l not in ("--fasta-search", "--predictor")]
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=cfg)
            w.generate(*HIGH)
            bad = w.run_job(w.job, FAKE_REPORT_RUNS="sample_1")
            self.assertNotEqual(bad.returncode, 0, bad.stdout + bad.stderr)
            self.assertIn("sample_0", bad.stdout + bad.stderr)
            good = w.run_job(w.job, FAKE_REPORT_RUNS="sample_0,sample_1")
            self.assertEqual(good.returncode, 0, good.stdout + good.stderr)


class StaleArtefactTests(unittest.TestCase):
    """A guard that checks a FILE cannot tell this search's output from the previous one's.

    Re-running a search into the same --out is routine (SKILL.md triggers: "re-run this
    search", "re-search with different parameters"). If DIA-NN then exits 0 without writing --
    no .NET on the node, an unmounted path -- the previous report is still there, the guards
    read it, and the job reports success on results from other parameters. Reproduced on HIVE
    with DIA-NN 2.7.0 (review srun 23512013): an old report.parquet + report.stats.tsv, a search
    with no DOTNET_ROOT, "ERROR: cannot read .raw files", exit 0, both files unchanged by md5,
    and check_report_runs.py printed "OK: report holds all 2 runs"."""

    def test_an_old_report_does_not_pass_for_a_search_that_wrote_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*HIGH)
            first = w.run_job(w.search_job, FAKE_REPORT_RUNS="sample_0,sample_1")
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            again = w.run_job(w.search_job)               # DIA-NN "succeeds", writes nothing
            self.assertNotEqual(again.returncode, 0, "passed on the previous search's report:\n"
                                + again.stdout + again.stderr)
            self.assertIn("report.parquet", again.stderr)

    def test_an_old_stats_file_does_not_stand_in_for_this_search(self):
        """The checker's no-parquet-reader route reads <report>.stats.tsv, so a stale one is the
        same hole. pyarrow and polars are hidden from the job to force that route."""
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*HIGH)
            hide = os.path.join(d, "no_parquet_readers")
            for mod in ("pyarrow", "polars"):
                os.makedirs(os.path.join(hide, mod))
                _write(os.path.join(hide, mod, "__init__.py"), "raise ImportError('hidden')\n")
            env = {"PYTHONPATH": hide}
            first = w.run_job(w.search_job, FAKE_REPORT_RUNS="sample_0,sample_1", **env)
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            self.assertIn("stats", first.stdout)
            # this search loses sample_1 and writes no stats file (--no-stats)
            again = w.run_job(w.search_job, FAKE_REPORT_RUNS="sample_0", FAKE_NO_STATS=1, **env)
            self.assertNotEqual(again.returncode, 0, "the previous search's stats file vouched "
                                "for this one:\n" + again.stdout + again.stderr)

    def test_an_old_library_does_not_pass_for_a_library_job_that_wrote_nothing(self):
        """Otherwise afterok releases the search against the PREVIOUS search's library."""
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*HIGH)
            self.assertEqual(w.run_job(w.lib_job, FAKE_WRITE_LIB=1).returncode, 0)
            again = w.run_job(w.lib_job, FAKE_WRITE_LIB=0)
            self.assertNotEqual(again.returncode, 0, again.stdout + again.stderr)
            self.assertIn("diann_lib.predicted.speclib", again.stderr)

    def test_the_single_job_path_does_not_pass_on_an_old_report(self):
        cfg = [l for l in LIBFREE_CFG if l not in ("--fasta-search", "--predictor")]
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=cfg)
            w.generate(*HIGH)
            first = w.run_job(w.job, FAKE_REPORT_RUNS="sample_0,sample_1")
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            again = w.run_job(w.job)
            self.assertNotEqual(again.returncode, 0, again.stdout + again.stderr)

    def test_an_inline_search_does_not_pass_on_an_old_report(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            first = w.run_inline(FAKE_WRITE_LIB=1, FAKE_REPORT_RUNS="sample_0,sample_1")
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            again = w.run_inline(FAKE_WRITE_LIB=1)
            self.assertNotEqual(again.returncode, 0, again.stdout + again.stderr)
            self.assertIn("report.parquet", again.stderr)

    def test_an_inline_search_does_not_pass_on_an_old_library(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            first = w.run_inline(FAKE_WRITE_LIB=1, FAKE_REPORT_RUNS="sample_0,sample_1")
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            again = w.run_inline(FAKE_WRITE_LIB=0, FAKE_REPORT_RUNS="sample_0,sample_1")
            self.assertNotEqual(again.returncode, 0, again.stdout + again.stderr)
            self.assertIn("predicted", again.stderr)


class DuplicateRunNameTests(unittest.TestCase):
    def test_inputs_sharing_a_run_name_are_refused_before_any_job_is_written(self):
        """DIA-NN names a run by its file name without the folder, so /a/s1.raw and /b/s1.raw
        collide. That is knowable from the input list, so it must stop generation -- not fail
        the search after it has used its SLURM hours."""
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            for sub in ("plate1", "plate2"):
                os.makedirs(os.path.join(d, sub))
            w.runs = [_write(os.path.join(d, sub, "sample_0.mzML"), "")
                      for sub in ("plate1", "plate2")]
            p = w.generate(*HIGH, check=False)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("sample_0", p.stderr)
            self.assertFalse(os.path.exists(w.lib_job))
            self.assertFalse(os.path.exists(w.search_job))


class ReportCheckTests(unittest.TestCase):
    """check_report_runs.py: the one place that decides whether a report holds every run."""

    def _stats(self, d, rows, name="report"):
        p = os.path.join(d, f"{name}.stats.tsv")
        with open(p, "w") as fh:
            fh.write("File.Name\tPrecursors.Identified\tProteins.Identified\n")
            for f, n in rows:
                fh.write(f"{f}\t{n}\t{n // 10}\n")
        return p

    def test_the_stats_file_route_needs_no_pyarrow(self):
        import check_report_runs as crr
        with tempfile.TemporaryDirectory() as d:
            files = [os.path.join(d, "a.raw"), os.path.join(d, "b.raw")]
            report = _write(os.path.join(d, "report.parquet"), "x")
            self._stats(d, [(files[0], 900), (files[1], 700)])
            ok, msg = crr.verify(report, files, parquet_runs=lambda p: None)
            self.assertTrue(ok, msg)
            self.assertIn("stats", msg)

    def test_a_zero_id_row_in_the_stats_file_counts_as_missing(self):
        """DIA-NN's stats file lists a run it failed on, with nothing identified; its report
        has no rows for it. Present-in-stats is not present-in-report."""
        import check_report_runs as crr
        with tempfile.TemporaryDirectory() as d:
            files = [os.path.join(d, "a.raw"), os.path.join(d, "b.raw")]
            report = _write(os.path.join(d, "report.parquet"), "x")
            self._stats(d, [(files[0], 900), (files[1], 0)])
            ok, msg = crr.verify(report, files, parquet_runs=lambda p: None)
            self.assertFalse(ok)
            self.assertIn("b", msg)

    def test_it_fails_closed_when_nothing_can_be_read(self):
        """No parquet reader and no stats file (--no-stats) means the run count is unknown.
        Passing that as success is the defect this check exists to remove."""
        import check_report_runs as crr
        with tempfile.TemporaryDirectory() as d:
            report = _write(os.path.join(d, "report.parquet"), "x")
            ok, msg = crr.verify(report, [os.path.join(d, "a.d")], parquet_runs=lambda p: None)
            self.assertFalse(ok)
            self.assertIn("pyarrow", msg)

    def test_inputs_sharing_a_run_name_fail(self):
        """DIA-NN's Run column is the file name without its path, so two inputs named alike in
        different folders become one Run: two samples silently merged for the DE step."""
        import check_report_runs as crr
        with tempfile.TemporaryDirectory() as d:
            files = [os.path.join(d, "x", "s1.raw"), os.path.join(d, "y", "s1.raw")]
            report = _write(os.path.join(d, "report.parquet"), "x")
            ok, msg = crr.verify(report, files, parquet_runs=lambda p: {"s1"})
            self.assertFalse(ok)
            self.assertIn("s1", msg)

    def test_a_missing_run_is_not_blamed_on_a_load_failure_without_evidence(self):
        """A run can be absent from the report because DIA-NN could not load it OR because it
        loaded and nothing passed the q-value filter (a blank, a failed injection). With no
        stats file to tell them apart, the message must not pick one."""
        import check_report_runs as crr
        with tempfile.TemporaryDirectory() as d:
            report = _write(os.path.join(d, "report.parquet"), "x")
            files = [os.path.join(d, "s1.raw"), os.path.join(d, "blank_01.raw")]
            ok, msg = crr.verify(report, files, parquet_runs=lambda p: {"s1"})
            self.assertFalse(ok)
            self.assertIn("blank_01", msg)
            self.assertIn("could not load", msg)
            self.assertIn("identified nothing", msg)

    def _signal_stats(self, d, rows):
        p = os.path.join(d, "report.stats.tsv")
        with open(p, "w") as fh:
            fh.write("File.Name\tPrecursors.Identified\tProteins.Identified\tTotal.Quantity"
                     "\tMS1.Signal\tMS2.Signal\n")
            for f, n, ms1, ms2 in rows:
                fh.write(f"{f}\t{n}\t{n // 10}\t{n * 1000}\t{ms1}\t{ms2}\n")

    def test_a_run_with_signal_but_no_identifications_is_reported_as_read(self):
        """DIA-NN's stats row separates the two when there is signal: a run it read has
        MS1/MS2 signal even with nothing identified. The job still fails -- a run with no
        identifications is not a sample the DE step can use -- but it says why."""
        import check_report_runs as crr
        with tempfile.TemporaryDirectory() as d:
            report = _write(os.path.join(d, "report.parquet"), "x")
            files = [os.path.join(d, "s1.raw"), os.path.join(d, "blank_01.raw"),
                     os.path.join(d, "broken.raw")]
            self._signal_stats(d, [(files[0], 2400, 2.2e12, 1.3e12),
                                   (files[1], 0, 3.1e10, 8.4e9),
                                   (files[2], 0, 0, 0)])
            ok, msg = crr.verify(report, files, parquet_runs=lambda p: {"s1"})
            self.assertFalse(ok)
            line = next(l for l in msg.splitlines() if l.strip().startswith("blank_01:"))
            self.assertIn("read", line)
            self.assertNotIn("could not load", line)
            line = next(l for l in msg.splitlines() if l.strip().startswith("broken:"))
            self.assertIn("could not load", line)

    def test_names_are_the_file_name_without_path_or_extension(self):
        import check_report_runs as crr
        self.assertEqual(crr.run_name("/data/a/Ex_42.raw"), "Ex_42")
        self.assertEqual(crr.run_name("/data/b/tims_1.d/"), "tims_1")
        self.assertEqual(crr.run_name("C:/x/y.mzML"), "y")

    @unittest.skipUnless(HAVE_PYARROW, "pyarrow not installed (CI is stdlib-only)")
    def test_the_parquet_route_reads_distinct_runs(self):
        import pyarrow as pa
        import pyarrow.parquet as pq
        import check_report_runs as crr
        with tempfile.TemporaryDirectory() as d:
            report = os.path.join(d, "report.parquet")
            pq.write_table(pa.table({"Run": ["a", "a", "b"], "Protein.Group": ["P"] * 3}), report)
            self.assertEqual(crr.runs_from_parquet(report), {"a", "b"})
            ok, _ = crr.verify(report, [os.path.join(d, "a.raw"), os.path.join(d, "b.raw")])
            self.assertTrue(ok)
            bad, msg = crr.verify(report, [os.path.join(d, n) for n in ("a.raw", "b.raw", "c.raw")])
            self.assertFalse(bad)
            self.assertIn("c", msg)


class RtProfilingTests(unittest.TestCase):
    """--reanalyse builds an empirical library in the first MBR pass, and --rt-profiling sets
    how that library is built. The search job is where that happens, so it keeps the flag."""

    def _search_tokens(self, cfg_lines):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=cfg_lines)
            w.generate(*HIGH)
            return _read(w.search_job).split()

    def test_the_search_job_keeps_rt_profiling_from_the_cfg(self):
        toks = self._search_tokens(LIBFREE_CFG)
        self.assertEqual(toks.count("--rt-profiling"), 1)
        self.assertEqual(toks.count("--reanalyse"), 1)

    def test_it_is_not_invented_when_the_cfg_lacks_it(self):
        toks = self._search_tokens([l for l in LIBFREE_CFG if l != "--rt-profiling"])
        self.assertNotIn("--rt-profiling", toks)


class ChainResourceTests(unittest.TestCase):
    """Steps 3 and 5 default to 64 CPUs, which rarely schedule on a congested preemptible
    queue; run_search.py gave no way to ask for less."""

    def test_chain_cpu_and_time_limits_pass_through(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], n_runs=3, ext=".d")
            a = argparse.Namespace(partition="low", account="publicgrp", qos="publicgrp-low-qos",
                                   max_simultaneous=None, assembly_cpus=24, libpred_cpus=8,
                                   time_per_file=5, assembly_mem=48)
            os.makedirs(w.out, exist_ok=True)
            run_search.run_diann_parallel(w.diann, w.cfg, w.runs, w.fasta, w.out, 16, a)
            s1 = _headers(os.path.join(w.out, "step1_libpred.sbatch"))
            s2 = _headers(os.path.join(w.out, "step2_firstpass.sbatch"))
            s3 = _headers(os.path.join(w.out, "step3_assembly.sbatch"))
            s5 = _headers(os.path.join(w.out, "step5_report.sbatch"))
            self.assertIn("#SBATCH --cpus-per-task=8", s1)
            self.assertIn("#SBATCH --time=5:00:00", s2)
            self.assertIn("#SBATCH --cpus-per-task=24", s3)
            self.assertIn("#SBATCH --cpus-per-task=24", s5)
            # fewer CPUs on a busy queue is not enough if the job still asks for 128 GB
            self.assertIn("#SBATCH --mem=48G", s3)
            self.assertIn("#SBATCH --mem=48G", s5)

    def test_run_search_exposes_the_flags(self):
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"), "--help"],
                           capture_output=True, text=True)
        for flag in ("--assembly-cpus", "--assembly-mem", "--libpred-cpus", "--time-per-file"):
            self.assertIn(flag, p.stdout)

    def test_chain_sizing_on_a_single_shot_search_is_not_silently_ignored(self):
        """Routing picks single-shot at <= 5 files, where these flags size nothing."""
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            p = w.generate(*HIGH, "--assembly-cpus", "8")
            self.assertIn("--assembly-cpus", p.stderr)
            self.assertIn("single-shot", p.stderr)

    def test_unset_flags_leave_the_generator_defaults_alone(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], n_runs=3, ext=".d")
            a = argparse.Namespace(partition="high", account="genome-center-grp", qos=None,
                                   max_simultaneous=None)
            os.makedirs(w.out, exist_ok=True)
            run_search.run_diann_parallel(w.diann, w.cfg, w.runs, w.fasta, w.out, 16, a)
            self.assertIn("#SBATCH --cpus-per-task=64",
                          _headers(os.path.join(w.out, "step3_assembly.sbatch")))


class ChainStaleArtefactTests(unittest.TestCase):
    """The 5-step chain's guards had the same blind spot as the single-shot jobs: must_exist()
    and step 5's `ls quant_step4/*.quant | wc -l` count files, not files THIS search wrote.
    Re-run a chain into the same --out (or resubmit a step, as references/watcher.md says to)
    and a step whose DIA-NN exits 0 having written nothing passes on the previous run's
    library, .quant or report."""

    def _chain(self, d):
        w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], n_runs=3, ext=".d")
        a = argparse.Namespace(partition="high", account="genome-center-grp", qos=None,
                               max_simultaneous=None)
        os.makedirs(w.out, exist_ok=True)
        run_search.run_diann_parallel(w.diann, w.cfg, w.runs, w.fasta, w.out, 16, a)
        return w

    def _old(self, w, *rel):
        for r in rel:
            p = os.path.join(w.out, r)
            os.makedirs(os.path.dirname(p), exist_ok=True)
            _write(p, "from the previous search")

    def _step(self, w, name, **env):
        return w.run_job(os.path.join(w.out, name), **env)

    def test_step1_does_not_pass_on_an_old_predicted_library(self):
        with tempfile.TemporaryDirectory() as d:
            w = self._chain(d)
            self._old(w, "step1.predicted.speclib")
            p = self._step(w, "step1_libpred.sbatch")
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)

    def test_step2_does_not_pass_on_an_old_quant(self):
        with tempfile.TemporaryDirectory() as d:
            w = self._chain(d)
            self._old(w, "quant_step2/sample_0.quant")
            p = self._step(w, "step2_firstpass.sbatch", SLURM_ARRAY_TASK_ID=0)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)

    def test_step3_does_not_pass_on_an_old_empirical_library(self):
        with tempfile.TemporaryDirectory() as d:
            w = self._chain(d)
            self._old(w, "empirical.parquet")
            p = self._step(w, "step3_assembly.sbatch")
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)

    def test_step4_does_not_pass_on_an_old_quant(self):
        with tempfile.TemporaryDirectory() as d:
            w = self._chain(d)
            self._old(w, "empirical.parquet", "quant_step2/sample_0.quant",
                      "quant_step4/sample_0.quant")
            p = self._step(w, "step4_finalpass.sbatch", SLURM_ARRAY_TASK_ID=0)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)

    def test_step5_does_not_pass_on_an_old_report(self):
        with tempfile.TemporaryDirectory() as d:
            w = self._chain(d)
            self._old(w, "report.parquet", "report.stats.tsv",
                      *[f"quant_step4/sample_{i}.quant" for i in range(3)])
            p = self._step(w, "step5_report.sbatch")               # DIA-NN writes nothing
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)

    def test_step5_counts_this_searchs_runs_not_every_quant_in_the_folder(self):
        """sample_2's step-4 task left no .quant, and a previous search's other_run.quant
        makes the folder count 3 of 3 anyway."""
        with tempfile.TemporaryDirectory() as d:
            w = self._chain(d)
            self._old(w, "quant_step4/sample_0.quant", "quant_step4/sample_1.quant",
                      "quant_step4/other_run.quant")
            p = self._step(w, "step5_report.sbatch",
                           FAKE_REPORT_RUNS="sample_0,sample_1,sample_2")
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("sample_2", p.stdout + p.stderr)

    def test_step5_passes_when_every_run_has_its_quant(self):
        with tempfile.TemporaryDirectory() as d:
            w = self._chain(d)
            self._old(w, *[f"quant_step4/sample_{i}.quant" for i in range(3)])
            p = self._step(w, "step5_report.sbatch",
                           FAKE_REPORT_RUNS="sample_0,sample_1,sample_2")
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("all 3 runs", p.stdout)


class ChainRunNameCollisionTests(unittest.TestCase):
    """DIA-NN names a run -- and diann_parallel names its .quant -- by the file name alone, so
    /plate1/s1.d and /plate2/s1.d are ONE Run and TWO array tasks writing one .quant.

    run_search.main() refuses that, but diann_parallel.py is ALSO run directly (its own
    docstring names that as the way to size the chain by hand), and that route had no check at
    all. Step 5's backstop did not catch it either: the count was derived per INPUT LINE of
    file_list.txt, so the single surviving s1.quant was counted once per input and the chain
    reported "built from all N runs" for a report holding N-1."""

    def _chain(self, d, **kw):
        w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], ext=".d", **kw)
        a = argparse.Namespace(partition="high", account="genome-center-grp", qos=None,
                               max_simultaneous=None)
        os.makedirs(w.out, exist_ok=True)
        run_search.run_diann_parallel(w.diann, w.cfg, w.runs, w.fasta, w.out, 16, a)
        return w

    def test_diann_parallel_run_directly_refuses_colliding_run_names(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], ext=".d")
            runs = []
            for sub in ("plate1", "plate2"):
                os.makedirs(os.path.join(d, sub))
                runs.append(_write(os.path.join(d, sub, "s1.d"), ""))
            out = os.path.join(d, "chain_out")
            p = subprocess.run(
                [sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                 "--diann", w.diann, "--raw", *runs, "--fasta", w.fasta,
                 "--out", out, "--cfg", w.cfg], capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("s1", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(out, "step4_finalpass.sbatch")),
                             "a chain was written for inputs that share a run name")

    def test_step5_counts_quant_files_not_input_lines(self):
        """Two inputs, one .quant they both map to. Counting input lines gives NQ=2 of 2 and
        prints "OK: report built from all 2 runs" over a report in which the two samples were
        merged -- exactly the silent dropped sample this guard exists to prevent."""
        with tempfile.TemporaryDirectory() as d:
            w = self._chain(d, n_runs=2)
            # a file_list.txt that collides: generated before main() refused it, or edited by
            # hand afterwards (references/watcher.md says to resubmit individual steps)
            _write(os.path.join(w.out, "file_list.txt"),
                   "/plate1/sample_0.d\n/plate2/sample_0.d\n")
            os.makedirs(os.path.join(w.out, "quant_step4"), exist_ok=True)
            _write(os.path.join(w.out, "quant_step4", "sample_0.quant"), "one file, two inputs")
            p = w.run_job(os.path.join(w.out, "step5_report.sbatch"),
                          FAKE_REPORT_RUNS="sample_0,sample_1")
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertNotIn("OK: report built from all", p.stdout)
            self.assertIn("distinct run names", p.stdout + p.stderr)

    def test_step5_still_passes_on_a_complete_chain(self):
        """The dedupe must not cost a legitimate chain its OK."""
        with tempfile.TemporaryDirectory() as d:
            w = self._chain(d, n_runs=3)
            os.makedirs(os.path.join(w.out, "quant_step4"), exist_ok=True)
            for i in range(3):
                _write(os.path.join(w.out, "quant_step4", f"sample_{i}.quant"), "q")
            p = w.run_job(os.path.join(w.out, "step5_report.sbatch"),
                          FAKE_REPORT_RUNS="sample_0,sample_1,sample_2")
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("all 3 runs", p.stdout)


class OutPathInjectionTests(unittest.TestCase):
    """clear_stale() and must_exist() put DOUBLE quotes around a path on purpose -- an array
    task's path carries $QUANT -- and bash re-reads `$` inside them. A `$(...)` in --out is
    therefore command substitution that runs when the JOB runs, inside the `rm -f --` the job
    performs before DIA-NN starts, so it also moves what gets deleted."""

    def test_run_search_refuses_an_out_that_would_execute(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            pwned = os.path.join(d, "pwned")
            w.out = os.path.join(d, f"o$(touch {pwned})")
            p = w.generate(*HIGH, check=False)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("REFUSING", p.stderr)
            self.assertFalse(os.path.exists(w.search_job), "a job was written anyway")
            self.assertFalse(os.path.exists(pwned))

    def test_the_chain_generator_refuses_it_too(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], n_runs=2, ext=".d")
            out = os.path.join(d, "o$(touch {})".format(os.path.join(d, "pwned")))
            p = subprocess.run(
                [sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                 "--diann", w.diann, "--raw", *w.runs, "--fasta", w.fasta,
                 "--out", out, "--cfg", w.cfg], capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("REFUSING", p.stderr)

    def test_only_the_characters_bash_re_reads_are_refused(self):
        """A path with spaces or parentheses is ordinary and stays inside the double quotes."""
        for ok in ("/data/run (2)/out", "/data/a b/out", "/data/o'ut"):
            self.assertEqual(diann_parallel.refuse_unsafe_path(ok), ok)
        for bad in ("/o$(id)", "/o`id`", '/o"x', "/o\\x", "/o\nx"):
            with self.assertRaises(SystemExit, msg=bad):
                diann_parallel.refuse_unsafe_path(bad)


class SingleShotStripTests(unittest.TestCase):
    """The flags the single-shot search job supplies ITSELF have to come out of the cfg, or
    they land on the command line twice."""

    def test_the_strip_list_is_the_chains_minus_a_stated_keep_list(self):
        """It was a hand-written list and had drifted from the chain's: --out, --f, --fasta,
        --threads, --lib and --out-lib were missing from it. Derived now, so a flag added to
        STRIP cannot go missing here in silence."""
        self.assertEqual(set(diann_parallel.STRIP) - set(run_search.SINGLE_SHOT_SEARCH_STRIP),
                         set(run_search.SINGLE_SHOT_SEARCH_KEEP))
        for f in ("--out", "--f", "--fasta", "--threads", "--lib", "--out-lib"):
            self.assertIn(f, run_search.SINGLE_SHOT_SEARCH_STRIP)
        # the job supplies none of these, so the cfg is their only source
        for f in ("--rt-profiling", "--xic", "--mobilograms", "--temp", "--no-norm"):
            self.assertNotIn(f, run_search.SINGLE_SHOT_SEARCH_STRIP)

    def test_a_cfg_that_names_its_own_out_does_not_reach_the_search_command(self):
        """The dangerous one: a second --out is a report DIA-NN may write somewhere the job
        neither clears beforehand nor checks afterwards, while the pre-delete still destroys
        the report at the path the job believes in."""
        with tempfile.TemporaryDirectory() as d:
            elsewhere = os.path.join(d, "elsewhere", "other.parquet")
            w = _Workspace(d, cfg_lines=LIBFREE_CFG + [
                f"--out {elsewhere}", "--threads 99",
                f"--lib {os.path.join(d, 'other.speclib')}"])
            w.generate(*HIGH)
            job = _read(w.search_job)
            line = next(l for l in job.splitlines() if l.startswith(w.diann + " "))
            toks = line.split()
            self.assertEqual(toks.count("--out"), 1, line)
            self.assertEqual(toks[toks.index("--out") + 1],
                             shlex.quote(os.path.join(w.out, "report.parquet")))
            self.assertEqual(toks.count("--threads"), 1, line)
            self.assertEqual(toks[toks.index("--threads") + 1], "4")
            self.assertEqual(toks.count("--lib"), 1, line)
            self.assertNotIn(elsewhere, job)

    def test_the_flags_the_job_does_not_supply_still_come_through(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=LIBFREE_CFG + [f"--temp {os.path.join(d, 'tmp')}"])
            w.generate(*HIGH)
            toks = _read(w.search_job).split()
            for f in ("--rt-profiling", "--xic", "--mobilograms", "--temp"):
                self.assertEqual(toks.count(f), 1, f)


class EmptyReportTests(unittest.TestCase):
    """DIA-NN creating report.parquet and then dying leaves it at 0 bytes. The sbatch route
    tests `[ -s ]`; the inline route tested os.path.exists, so the empty file passed and
    reached the parquet reader."""

    def test_the_sbatch_route_fails_on_a_zero_byte_report(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*HIGH)
            p = w.run_job(w.search_job, FAKE_EMPTY_REPORT=1)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)

    def test_the_inline_route_fails_on_it_the_same_way(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            p = w.run_inline(FAKE_WRITE_LIB=1, FAKE_EMPTY_REPORT=1)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("missing or empty", p.stdout + p.stderr)
            self.assertEqual(os.path.getsize(os.path.join(w.out, "report.parquet")), 0)


class InlineSetAsideTests(unittest.TestCase):
    """The sbatch route can only DELETE the previous artefacts (clear_stale runs on a compute
    node long after generation). Inline we are the ones re-running, so they are renamed the way
    --sbatch already renames a job script it is about to replace (`existing_file_moved_to`):
    a re-run whose DIA-NN dies then still leaves the user the report they had."""

    def test_a_re_run_sets_the_previous_report_aside_instead_of_deleting_it(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            first = w.run_inline(FAKE_WRITE_LIB=1, FAKE_REPORT_RUNS="sample_0,sample_1")
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            report = os.path.join(w.out, "report.parquet")
            with open(report, "rb") as fh:
                kept = fh.read()
            again = w.run_inline(FAKE_WRITE_LIB=1)        # "succeeds", writes nothing
            self.assertNotEqual(again.returncode, 0, again.stdout + again.stderr)
            stale = [f for f in os.listdir(w.out) if f.startswith("report.parquet.stale-")]
            self.assertEqual(len(stale), 1, os.listdir(w.out))
            with open(os.path.join(w.out, stale[0]), "rb") as fh:
                self.assertEqual(fh.read(), kept)
            self.assertFalse(os.path.exists(report), "the stale report was left in place")


class InputListingTests(unittest.TestCase):
    """report_guard() bakes the input list into the job script, so the list has to be versioned
    with it. One shared search_input_files.txt meant generating a second search into the same
    --out rewrote the list the FIRST, still unsubmitted, job would read: that job would then
    check its report against a cohort it never searched."""

    def _generate(self, w, job, runs):
        return subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "run_search.py"),
             "--tools", w.tools, "--bundle", w.bundle, "--params", w.cfg,
             "--fasta", w.fasta, "--out", w.out, "--files", *runs,
             "--threads", "4", "--sbatch", job, *HIGH], capture_output=True, text=True)

    def test_the_input_list_is_named_after_the_job_that_reads_it(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, n_runs=3)
            p = self._generate(w, os.path.join(d, "job_a.sh"), w.runs[:2])
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            list_a = os.path.join(w.out, "job_a_input_files.txt")
            self.assertIn(list_a, _read(os.path.join(d, "job_a_2_search.sh")))
            self.assertEqual(_read(list_a).split(), w.runs[:2])

            p = self._generate(w, os.path.join(d, "job_b.sh"), w.runs)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(_read(list_a).split(), w.runs[:2],
                             "generating a second job rewrote the first job's input list")
            self.assertEqual(_read(os.path.join(w.out, "job_b_input_files.txt")).split(),
                             w.runs)


class _FakeSacctmgr:
    """Put a sacctmgr on PATH that prints `assoc`, for the duration of `test`."""

    def __init__(self, test, assoc, user="tester"):
        tmp = tempfile.TemporaryDirectory()
        test.addCleanup(tmp.cleanup)
        binp = os.path.join(tmp.name, "bin")
        os.makedirs(binp)
        _write(os.path.join(binp, "sacctmgr"),
               "#!/bin/sh\ncat <<'ASSOC'\n" + assoc + "ASSOC\n", 0o755)
        env = mock.patch.dict(os.environ, {"PATH": binp + os.pathsep + os.environ["PATH"],
                                           "USER": user})
        env.start()
        test.addCleanup(env.stop)


class AmbiguousAssociationTests(unittest.TestCase):
    """`(preferred or hits)[0]` picked whichever association sacctmgr happened to list first,
    with no output at all. For a user in two labs that decides who is BILLED for a multi-hour
    search."""

    TWO_LABS = "labA|low|labA-low-qos\nlabB|low|labB-low-qos\n"
    ONE_LAB = "labA|low|labA-low-qos\nlabA|gpu|labA-gpu-qos\n"

    def test_a_partition_matching_two_accounts_is_refused(self):
        _FakeSacctmgr(self, self.TWO_LABS)
        with self.assertRaises(SystemExit) as e:
            run_search.slurm_queue(partition="low")
        msg = str(e.exception)
        self.assertIn("labA", msg)
        self.assertIn("labB", msg)
        self.assertIn("billed", msg)

    def test_one_account_on_two_partitions_says_which_it_took(self):
        """Same account either way, so this is a scheduling choice, not a billing one -- but
        it is still a choice the user did not make, so it is not made in silence."""
        _FakeSacctmgr(self, self.ONE_LAB)
        err = io.StringIO()
        with mock.patch("sys.stderr", err):
            part, acct, qos = run_search.slurm_queue(account="labA")
        self.assertEqual(acct, "labA")
        self.assertIn("WARNING", err.getvalue())
        self.assertIn("labA|low", err.getvalue())
        self.assertIn("labA|gpu", err.getvalue())

    def test_an_unambiguous_match_says_nothing(self):
        _FakeSacctmgr(self, HIVE_ASSOCIATIONS)
        err = io.StringIO()
        with mock.patch("sys.stderr", err):
            got = run_search.slurm_queue(partition="low")
        self.assertEqual(got, ("low", "publicgrp", "publicgrp-low-qos"))
        self.assertEqual(err.getvalue(), "")


class BlankAssociationFieldTests(unittest.TestCase):
    """A blank partition or QOS field means "no restriction". That is a fine DEFAULT, but it
    must not confirm a value the USER typed -- SLURM will reject it, hours later."""

    BLANK_QOS = "labA|low|\n"
    BLANK_PARTITION = "labA||labA-qos\n"

    def test_a_bogus_qos_is_not_confirmed_by_a_blank_qos_field(self):
        _FakeSacctmgr(self, self.BLANK_QOS)
        with self.assertRaises(SystemExit) as e:
            run_search.slurm_queue(qos="totally-bogus-qos")
        self.assertIn("totally-bogus-qos", str(e.exception))

    def test_a_bogus_partition_is_not_confirmed_by_a_blank_partition_field(self):
        _FakeSacctmgr(self, self.BLANK_PARTITION)
        with self.assertRaises(SystemExit) as e:
            run_search.slurm_queue(partition="no-such-partition")
        self.assertIn("no-such-partition", str(e.exception))

    def test_a_blank_field_is_still_used_when_nothing_was_given(self):
        """Detection on a cluster laid out unlike HIVE must keep working."""
        _FakeSacctmgr(self, self.BLANK_PARTITION)
        self.assertEqual(run_search.slurm_queue(), (None, "labA", "labA-qos"))


class CompleteQueueValidationTests(unittest.TestCase):
    """A complete --partition/--account pair used to return BEFORE the associations were read,
    so the only overrides that got checked were the incomplete ones. It is still honoured --
    nothing has to be invented, and a reservation may exist sacctmgr does not show -- but a
    mismatch is said out loud instead of surfacing as a SLURM rejection at submit time."""

    def test_a_pair_no_association_has_is_honoured_but_reported(self):
        _FakeSacctmgr(self, HIVE_ASSOCIATIONS)
        err = io.StringIO()
        with mock.patch("sys.stderr", err):
            got = run_search.slurm_queue(partition="low", account="genome-center-grp")
        self.assertEqual(got, ("low", "genome-center-grp", None))
        self.assertIn("WARNING", err.getvalue())
        self.assertIn("publicgrp|low", err.getvalue())      # names what the user CAN use

    def test_a_pair_an_association_has_says_nothing(self):
        _FakeSacctmgr(self, HIVE_ASSOCIATIONS)
        err = io.StringIO()
        with mock.patch("sys.stderr", err):
            got = run_search.slurm_queue(partition="high", account="genome-center-grp")
        self.assertEqual(got, ("high", "genome-center-grp", None))
        self.assertEqual(err.getvalue(), "")


class ExecutableTests(unittest.TestCase):
    def test_the_report_checker_is_executable(self):
        """CI requires every shebang script to ship 100755."""
        p = os.path.join(SCRIPTS, "check_report_runs.py")
        self.assertTrue(os.path.exists(p))
        self.assertTrue(os.stat(p).st_mode & stat.S_IXUSR)


if __name__ == "__main__":
    unittest.main(verbosity=2)
