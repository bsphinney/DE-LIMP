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
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
import unittest

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

    def generate(self, *queue):
        p = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "run_search.py"),
             "--tools", self.tools, "--bundle", self.bundle, "--params", self.cfg,
             "--fasta", self.fasta, "--out", self.out, "--files", *self.runs,
             "--threads", "4", "--sbatch", self.job, *queue],
            capture_output=True, text=True)
        if p.returncode != 0:
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
                                   time_per_file=5)
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

    def test_run_search_exposes_the_flags(self):
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"), "--help"],
                           capture_output=True, text=True)
        for flag in ("--assembly-cpus", "--libpred-cpus", "--time-per-file"):
            self.assertIn(flag, p.stdout)

    def test_unset_flags_leave_the_generator_defaults_alone(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], n_runs=3, ext=".d")
            a = argparse.Namespace(partition="high", account="genome-center-grp", qos=None,
                                   max_simultaneous=None)
            os.makedirs(w.out, exist_ok=True)
            run_search.run_diann_parallel(w.diann, w.cfg, w.runs, w.fasta, w.out, 16, a)
            self.assertIn("#SBATCH --cpus-per-task=64",
                          _headers(os.path.join(w.out, "step3_assembly.sbatch")))


class ExecutableTests(unittest.TestCase):
    def test_the_report_checker_is_executable(self):
        """CI requires every shebang script to ship 100755."""
        p = os.path.join(SCRIPTS, "check_report_runs.py")
        self.assertTrue(os.path.exists(p))
        self.assertTrue(os.stat(p).st_mode & stat.S_IXUSR)


if __name__ == "__main__":
    unittest.main(verbosity=2)
