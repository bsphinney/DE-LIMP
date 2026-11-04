#!/usr/bin/env python3
"""
A DIA-NN --dda search whose SECOND pass (match-between-runs) crashes falls back to the complete
first-pass report DIA-NN already wrote, and says so (I50).

A Core submission (2026-10-29, DIA-NN 2.7.0 --dda, two timsTOF .d): the second pass crashed, and the
search folder held `report-first-pass.parquet` (with its .manifest.txt, .pr_matrix.tsv and
.stats.tsv) and no report.parquet, report.pg_matrix.tsv or report.stats.tsv. The job failed, the
skill had nothing to analyse, and the session left the skill's path. Now the search job copies the
first pass onto the report (dda_first_pass.record), records it (dda_first_pass.json and
search_provenance.json `report_used: "first_pass"`), and the caveat -- ONE definition,
dda_first_pass.caution() -- reaches the Methods, AUDIT.md and the watcher.

What must not change: a DIA search's job (no fallback code at all); a crash with no first pass
(fails with DIA-NN's own exit status); a final report that exists (never replaced, whatever the
exit status). A search provenance from before 2.14.2 has no `report_used`, and the watcher shows
no caution for it (Review Focus 1). A path with a space or a quote stays one word (Review Focus 2).

Synthetic files only, with the layout read on HIVE; fakes stand in for DIA-NN and SLURM.
"""
import hashlib
import json
import os
import shlex
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
WATCH = os.path.join(SCRIPTS, "watch_run.sh")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

from job_env import job_env  # noqa: E402  (env for running generated scripts)

import dda_first_pass  # noqa: E402
import audit_results as ar  # noqa: E402
import make_methods as mm  # noqa: E402
# Module imports, not `from ... import`: a TestCase bound here would be collected twice.
import test_step1b_mass_accuracy_probe as fx  # noqa: E402  (fake-executable helper)

# A DIA-NN stand-in. FAKE_MODE says what its SEARCH does (a call with --out-lib predicts the
# library and succeeds):
#   first_pass_132  writes <stem>-first-pass.parquet (+ .stats.tsv, .pr_matrix.tsv) and dies, 132
#   first_pass_0    the same, but exits 0 (DIA-NN exits 0 on some fatal errors)
#   none_132        writes nothing, exits 132
#   final_132       writes the final report AND the first pass, then exits 132
#   final_0         a normal search: the final report, exit 0
FAKE_DIANN = r"""#!/bin/bash
out=""; outlib=""; prev=""; files=()
for a in "$@"; do
  [ "$prev" = --f ] && files+=("$a")
  [ "$prev" = --out ] && out="$a"
  [ "$prev" = --out-lib ] && outlib="$a"
  prev="$a"
done
if [ -n "$outlib" ]; then echo lib > "${outlib%.*}.predicted.speclib"; exit 0; fi
stem="${out%.parquet}"
write() {   # $1 = report path to write, with its stats; the pass is in its rows
  python3 - "$1" "${files[@]}" <<'PY'
import os, sys
out, inputs = sys.argv[1], sys.argv[2:]
runs = [os.path.splitext(os.path.basename(f.rstrip("/")))[0] for f in inputs]
which = "first" if out.endswith("-first-pass.parquet") else "final"
try:
    import pyarrow as pa, pyarrow.parquet as pq
    pq.write_table(pa.table({"Run": runs, "Protein.Group": [which] * len(runs)}), out)
except ImportError:
    open(out, "wb").write(b"PAR1 placeholder: no pyarrow in this interpreter " + out.encode())
stem = out[:-len(".parquet")]
with open(stem + ".stats.tsv", "w") as fh:
    fh.write("File.Name\tPrecursors.Identified\tProteins.Identified\n")
    for f in inputs:
        fh.write(f"{f}\t1500\t150\n")
PY
}
case "${FAKE_MODE:-final_0}" in
  first_pass_132) write "$stem-first-pass.parquet"; echo pr > "$stem-first-pass.pr_matrix.tsv"
                  echo "Segmentation fault" >&2; exit 132 ;;
  first_pass_0)   write "$stem-first-pass.parquet"; exit 0 ;;
  none_132)       exit 132 ;;
  final_132)      write "$stem-first-pass.parquet"; write "$out"; exit 132 ;;
  *)              write "$out"; exit 0 ;;
esac
"""


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def run_job_snippet(d, mode, out_name="out", provenance=True):
    """The snippet job_lines() writes, run under the job's own `set -euo pipefail`, with the fake
    DIA-NN as the search. -> (rc, stdout+stderr, out dir, report, provenance)."""
    out = os.path.join(d, out_name)
    os.makedirs(out, exist_ok=True)
    fake = fx._exe(os.path.join(d, "diann"), FAKE_DIANN)
    report = os.path.join(out, "report.parquet")
    prov = os.path.join(out, "search_provenance.json")
    if provenance:
        with open(prov, "w") as fh:
            json.dump({"engine": "diann", "version": "2.7.0"}, fh, indent=2)
    raws = []
    for name in ("runA.d", "runB.d"):
        raws.append(os.path.join(d, name))
        os.makedirs(raws[-1], exist_ok=True)
    q = shlex.quote
    search_cmd = (f"{q(fake)} --dda " + " ".join(f"--f {q(r)}" for r in raws)
                  + f" --reanalyse --out {q(report)} --threads 4")
    body = "\n".join(["set -euo pipefail", f"cd {q(out)}",
                      *dda_first_pass.job_lines(search_cmd, report, prov),
                      'echo "REACHED THE GUARD"', ""])
    script = os.path.join(d, "job.sh")
    with open(script, "w") as fh:
        fh.write("#!/bin/bash\n" + body)
    r = subprocess.run(["bash", script], capture_output=True, text=True, timeout=120,
                       env=job_env(d, FAKE_MODE=mode))
    return r.returncode, r.stdout + r.stderr, out, report, prov


class FirstPassPath(unittest.TestCase):
    def test_the_name_dia_nn_2_7_0_writes(self):
        """A Core member's output/search_dda/: report-first-pass.parquet beside the report's stem."""
        self.assertEqual(dda_first_pass.first_pass_path("/s/out/report.parquet"),
                         "/s/out/report-first-pass.parquet")
        self.assertEqual(dda_first_pass.FIRST_PASS_SUFFIX, "-first-pass.parquet")


class JobSnippet(unittest.TestCase):
    def test_crash_with_a_first_pass_falls_back(self):
        with tempfile.TemporaryDirectory() as d:
            rc, log, out, report, provenance = run_job_snippet(d, "first_pass_132")
            self.assertEqual(rc, 0, log)
            self.assertIn("REACHED THE GUARD", log)
            self.assertIn("DIANN_DDA_FIRST_PASS_FALLBACK", log)
            self.assertTrue(os.path.getsize(report) > 0)
            fp = dda_first_pass.first_pass_path(report)
            self.assertEqual(sha256(report), sha256(fp))
            self.assertFalse(os.path.islink(report), "a copy, not a link: a link breaks in a zip")
            # the first pass's stats under the report's name: the guard's fallback evidence
            self.assertTrue(os.path.isfile(os.path.join(out, "report.stats.tsv")))
            # ...and its precursor matrix, which a delivery takes under the report's name
            # (core_submission.SEARCH_FILES), as it takes the report (L8R-F8)
            with open(os.path.join(out, "report.pr_matrix.tsv")) as fh:
                self.assertEqual(fh.read(), "pr\n")
            self.assertFalse(os.path.exists(report + ".part"))
            prov = json.load(open(provenance))
            self.assertEqual(prov["report_used"], "first_pass")
            self.assertEqual(prov["dda_first_pass"]["exit_code"], 132)
            self.assertEqual(prov["dda_first_pass"]["first_pass_sha256"], sha256(fp))
            self.assertEqual(prov["version"], "2.7.0", "the rest of the record is kept")
            rec = json.load(open(os.path.join(out, dda_first_pass.RECORD)))
            self.assertEqual(rec["report_used"], "first_pass")
            self.assertIn("at", rec)
            self.assertEqual(rec["caveat"], dda_first_pass.caution(prov))

    def test_dia_nn_exits_0_without_the_report_still_falls_back(self):
        """DIA-NN exits 0 on some fatal errors: the snippet keys on the missing report."""
        with tempfile.TemporaryDirectory() as d:
            rc, log, _, report, provenance = run_job_snippet(d, "first_pass_0")
            self.assertEqual(rc, 0, log)
            self.assertTrue(os.path.getsize(report) > 0)
            self.assertEqual(json.load(open(provenance))["dda_first_pass"]["exit_code"], 0)

    def test_crash_without_a_first_pass_still_fails(self):
        with tempfile.TemporaryDirectory() as d:
            rc, log, out, report, provenance = run_job_snippet(d, "none_132")
            self.assertEqual(rc, 132, log)
            self.assertNotIn("REACHED THE GUARD", log)
            self.assertFalse(os.path.exists(report))
            self.assertNotIn("report_used", json.load(open(provenance)))
            self.assertFalse(os.path.exists(os.path.join(out, dda_first_pass.RECORD)))

    def test_a_final_report_is_never_replaced(self):
        """The final report exists: the second pass wrote it. Its exit status decides, as before
        (2.14.1 failed the job on it), and the first pass never replaces it."""
        with tempfile.TemporaryDirectory() as d:
            rc, log, out, report, provenance = run_job_snippet(d, "final_132")
            self.assertEqual(rc, 132, log)
            self.assertNotEqual(sha256(report), sha256(dda_first_pass.first_pass_path(report)))
            self.assertNotIn("report_used", json.load(open(provenance)))
            self.assertFalse(os.path.exists(os.path.join(out, dda_first_pass.RECORD)))

    def test_a_normal_search_is_untouched(self):
        with tempfile.TemporaryDirectory() as d:
            rc, log, out, report, provenance = run_job_snippet(d, "final_0")
            self.assertEqual(rc, 0, log)
            self.assertNotIn("DIANN_DDA_FIRST_PASS_FALLBACK", log)
            self.assertNotIn("report_used", json.load(open(provenance)))

    def test_quoted_paths_with_spaces(self):
        """Review Focus 2: a service folder with a space or a quote stays one word."""
        with tempfile.TemporaryDirectory() as d:
            rc, log, out, report, provenance = run_job_snippet(d, "first_pass_132",
                                                               out_name="Smith lab's out")
            self.assertEqual(rc, 0, log)
            self.assertTrue(os.path.getsize(report) > 0)
            self.assertEqual(json.load(open(provenance))["report_used"], "first_pass")

    def test_a_job_whose_provenance_is_missing_fails_and_copies_nothing(self):
        """L8R-F9: the Methods, AUDIT.md and the watcher read only search_provenance.json, so a
        fallback the job cannot record there is a failed job, not a silent first pass."""
        with tempfile.TemporaryDirectory() as d:
            rc, log, out, report, provenance = run_job_snippet(d, "first_pass_132",
                                                               provenance=False)
            self.assertNotEqual(rc, 0, log)
            self.assertIn("search_provenance.json", log)
            self.assertFalse(os.path.exists(report))
            self.assertFalse(os.path.exists(os.path.join(out, dda_first_pass.RECORD)))

    def test_the_record_refuses_when_the_final_report_exists(self):
        with tempfile.TemporaryDirectory() as d:
            report = os.path.join(d, "report.parquet")
            fp = dda_first_pass.first_pass_path(report)
            for p in (report, fp):
                with open(p, "wb") as fh:
                    fh.write(os.path.basename(p).encode())
            with self.assertRaises(ValueError):
                dda_first_pass.record(report, fp, None, 132)
            self.assertEqual(open(report, "rb").read(), b"report.parquet")


class Caution(unittest.TestCase):
    def test_none_without_report_used(self):
        self.assertIsNone(dda_first_pass.caution({}))
        self.assertIsNone(dda_first_pass.caution(None))
        self.assertIsNone(dda_first_pass.caution({"report_used": "final"}))

    def test_filled_in(self):
        c = dda_first_pass.caution({"report_used": "first_pass",
                                    "dda_first_pass": {"exit_code": 132}})
        self.assertIn("without match-between-runs", c)
        self.assertIn("no identification was transferred between runs", c)

    def test_it_states_no_direction_nobody_measured_and_no_exit_code(self):
        """L8R-F4: there is no second-pass report to compare with, and on SET28 the second pass
        had FEWER precursors in 18 of 28 runs (pass_comparison.py), so "fewer precursors" and
        "runs agree less" were claims. A shell status is not Methods text either: it stays in
        dda_first_pass.json and AUDIT.md's evidence."""
        c = dda_first_pass.caution({"report_used": "first_pass",
                                    "dda_first_pass": {"exit_code": 132}})
        for claim in ("fewer", "agree less", "exit", "132"):
            self.assertNotIn(claim, c)
        self.assertIn("more values may be missing across runs", c)
        words = ar.PLAIN["dda_first_pass"][1]["WARN"]
        for claim in ("fewer", "agree less"):
            self.assertNotIn(claim, words)
        self.assertIn("found directly", words)


# A library-free cfg as estimate_params.py writes one, mass accuracy pinned (no probe)
LIBFREE = ("--qvalue 0.01", "--matrices", "--xic 10", "--mobilograms", "--fasta-search",
           "--gen-spec-lib", "--predictor", "--reanalyse", "--rt-profiling", "--cut K*,R*",
           "--mass-acc 20", "--mass-acc-ms1 10")


class Search:
    """A two-job DIA-NN search generated by run_search.py --sbatch, with FAKE_DIANN as DIA-NN:
    mzML inputs (no .NET, no multi-GB fixture on an HFS+ TMPDIR), raw data and --out in folders
    with a space, --out with a quote too (Review Focus 2)."""

    def __init__(self, d, dda=True, out_name="Smith lab's out"):
        self.d = d
        raw = os.path.join(d, "raw data")
        os.makedirs(raw)
        self.raws = []
        for name in ("runA.mzML", "runB.mzML"):
            self.raws.append(os.path.join(raw, name))
            with open(self.raws[-1], "w") as fh:
                fh.write("<mzML/>")
        self.fasta = os.path.join(d, "db.fasta")
        with open(self.fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDER\n")
        self.cfg = os.path.join(d, "diann.cfg")
        with open(self.cfg, "w") as fh:
            fh.write("\n".join(LIBFREE + (("--dda",) if dda else ())) + "\n")
        self.tools = os.path.join(d, "tools.json")
        with open(self.tools, "w") as fh:
            json.dump({"diann": fx._exe(os.path.join(d, "diann"), FAKE_DIANN)}, fh)
        self.bundle = os.path.join(d, "bundle.json")
        with open(self.bundle, "w") as fh:
            json.dump({"acquisition": "DDA" if dda else "DIA", "engine": {"name": "diann"}}, fh)
        self.out = os.path.join(d, out_name)

    def env(self, mode="final_0"):
        # no sbatch on PATH: not a login node (job_env has no SLURM at the fixed path either)
        return job_env(self.d, FAKE_MODE=mode, PATH=os.pathsep.join(
            [os.path.dirname(sys.executable), "/usr/bin", "/bin"]))

    def generate(self, *extra, sbatch=True):
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                            "--tools", self.tools, "--bundle", self.bundle, "--params", self.cfg,
                            "--fasta", self.fasta, "--out", self.out, "--files", *self.raws,
                            "--engine", "diann", "--threads", "4", *extra]
                           + (["--sbatch", os.path.join(self.d, "job.sh")] if sbatch else []),
                           capture_output=True, text=True, env=self.env(),
                           cwd=self.d, timeout=240)
        assert p.returncode == 0, p.stdout + p.stderr
        return p

    def run(self, script, mode):
        return subprocess.run(["bash", os.path.join(self.d, script)], capture_output=True,
                              text=True, env=self.env(mode), timeout=240)

    def prov(self):
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            return json.load(fh)


class TheSearchJob(unittest.TestCase):
    def test_dia_job_is_unchanged(self):
        """A DIA search's job carries none of the fallback (byte-identical to 2.14.1's)."""
        with tempfile.TemporaryDirectory() as d:
            w = Search(d, dda=False)
            w.generate()
            srch = open(os.path.join(d, "job_2_search.sh")).read()
            for marker in ("dda_first_pass", "_dia_rc", "set +e", "first-pass"):
                self.assertNotIn(marker, srch)
            lines = srch.splitlines()
            i = next(k for k, ln in enumerate(lines) if "--reanalyse" in ln)
            self.assertIn("the report", lines[i + 1], "the guard follows the search at once")

    def test_dda_job_falls_back_and_its_guard_passes(self):
        with tempfile.TemporaryDirectory() as d:
            w = Search(d)
            w.generate()
            srch = open(os.path.join(d, "job_2_search.sh")).read()
            self.assertEqual(len([ln for ln in srch.splitlines() if "--reanalyse" in ln]), 1)
            self.assertIn("dda_first_pass.py", srch)
            for script in ("job_1_lib.sh", "job_2_search.sh"):
                r = w.run(script, "first_pass_132")
                self.assertEqual(r.returncode, 0, script + "\n" + r.stdout + r.stderr)
            self.assertIn("DIANN_DDA_FIRST_PASS_FALLBACK", r.stdout + r.stderr)
            prov = w.prov()
            self.assertEqual(prov["report_used"], "first_pass")
            self.assertEqual(prov["engine"], "diann")
            self.assertTrue(os.path.getsize(os.path.join(w.out, "report.parquet")) > 0)

    def test_dda_job_without_a_first_pass_still_fails(self):
        with tempfile.TemporaryDirectory() as d:
            w = Search(d)
            w.generate()
            w.run("job_1_lib.sh", "none_132")
            r = w.run("job_2_search.sh", "none_132")
            self.assertEqual(r.returncode, 132, r.stdout + r.stderr)

    def test_a_resubmitted_job_that_succeeds_reads_as_a_second_pass(self):
        """L8R-F1: a fallback, then job_2_search.sh resubmitted as it is (what the two-job split
        is for: requeue against the same library) and this time the second pass succeeds. The
        provenance must no longer say first_pass -- Methods, AUDIT.md, the watcher and FRAN's
        manifest all read it."""
        with tempfile.TemporaryDirectory() as d:
            w = Search(d)
            w.generate()
            w.run("job_1_lib.sh", "final_0")
            r = w.run("job_2_search.sh", "first_pass_132")
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertEqual(w.prov()["report_used"], "first_pass")
            r = w.run("job_2_search.sh", "final_0")
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            prov = w.prov()
            self.assertNotIn("report_used", prov)
            self.assertNotIn("dda_first_pass", prov)
            self.assertIsNone(dda_first_pass.caution(prov))
            self.assertFalse(os.path.exists(os.path.join(w.out, dda_first_pass.RECORD)))
            self.assertEqual(prov["engine"], "diann", "every other key is kept")
            self.assertIn("result", prov)
            sp = os.path.join(w.out, "search_provenance.json")
            para = mm.search_paragraph(mm.search_record(None, sp, None))
            self.assertNotIn("first-pass report", para)

    def test_an_earlier_search_s_first_pass_files_are_cleared(self):
        """L8R-F8: every first-pass file of an earlier search goes before this one runs, its
        matrix and manifest too, so none of them can sit beside this search's report."""
        with tempfile.TemporaryDirectory() as d:
            w = Search(d)
            w.generate()
            os.makedirs(w.out, exist_ok=True)
            stale = [os.path.join(w.out, f) for f in (
                "report-first-pass.parquet", "report-first-pass.stats.tsv",
                "report-first-pass.pr_matrix.tsv", "report-first-pass.manifest.txt",
                "report.pr_matrix.tsv", dda_first_pass.RECORD)]
            for f in stale:
                with open(f, "w") as fh:
                    fh.write("an earlier search\n")
            w.run("job_1_lib.sh", "final_0")
            r = w.run("job_2_search.sh", "final_0")
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertEqual([f for f in stale if os.path.exists(f)], [])

    def test_inline_dda_crash_falls_back(self):
        """The inline route does the same, in Python, and main() records it."""
        with tempfile.TemporaryDirectory() as d:
            w = Search(d)
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                                "--tools", w.tools, "--bundle", w.bundle, "--params", w.cfg,
                                "--fasta", w.fasta, "--out", w.out, "--files", *w.raws,
                                "--engine", "diann", "--threads", "4"],
                               capture_output=True, text=True, env=w.env("first_pass_132"),
                               timeout=240)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            prov = w.prov()
            self.assertEqual(prov["report_used"], "first_pass")
            self.assertEqual(prov["dda_first_pass"]["exit_code"], 132)
            self.assertIsNotNone(dda_first_pass.caution(prov))


FALLBACK_PROV = {"engine": "diann", "version": "2.7.0", "search_mode": "single_shot",
                 "resolved_command": "/x/diann", "report_used": "first_pass",
                 "dda_first_pass": {"report_used": "first_pass", "exit_code": 132}}


class MethodsAndAudit(unittest.TestCase):
    def test_methods_and_audit_carry_the_caveat(self):
        caveat = dda_first_pass.caution(FALLBACK_PROV)
        with tempfile.TemporaryDirectory() as d:
            sp = os.path.join(d, "search_provenance.json")
            with open(sp, "w") as fh:
                json.dump(FALLBACK_PROV, fh)
            para = mm.search_paragraph(mm.search_record(None, sp, None))
            self.assertIn(caveat, para)
            self.assertNotIn("with match-between-runs.", para)
            # the Methods' parameter table says which report was quantified, never "MBR: on"
            row = mm.quantified_report_row(mm.search_record(None, sp, None), None)
            self.assertEqual(row[0], "Report quantified")
            self.assertIn("first-pass report", row[1])
            self.assertIn("report_used", row[2])
            self.assertNotIn("exit 132", para)
            findings = []
            ar.audit_dda_first_pass(findings, {"input": os.path.join(d, "report.parquet")})
            self.assertEqual([(f["check"], f["status"]) for f in findings],
                             [("dda_first_pass", "WARN")])
            self.assertTrue(findings[0]["message"].startswith(caveat), findings[0]["message"])
            # the shell status is the audit's evidence (AUDIT.md and its detail), not the Methods'
            self.assertIn("DIA-NN exited 132", findings[0]["message"])
            self.assertEqual(findings[0]["detail"]["exit_code"], 132)
            label, meaning = ar.client_row(findings[0])
            self.assertNotIn("_", label)
            self.assertNotEqual(meaning, ar.GENERIC_MEANING["WARN"])

    def test_an_ordinary_search_has_no_caveat(self):
        with tempfile.TemporaryDirectory() as d:
            sp = os.path.join(d, "search_provenance.json")
            with open(sp, "w") as fh:
                json.dump({k: v for k, v in FALLBACK_PROV.items()
                           if k not in ("report_used", "dda_first_pass")}, fh)
            para = mm.search_paragraph(mm.search_record(None, sp, None))
            self.assertNotIn("first-pass report", para)
            findings = []
            ar.audit_dda_first_pass(findings, {"input": os.path.join(d, "report.parquet")})
            self.assertEqual(findings, [])


class StatusViews(unittest.TestCase):
    """L8R-F3: the Core run registry (record_run.py), the job-end Slack post (notify_slack.py)
    and checkpoint.py status read the same caution (rule 3) -- and none says "Match-between-runs:
    on" for a search whose second pass failed."""

    def _out(self, d, prov=FALLBACK_PROV):
        out = os.path.join(d, "search out")
        os.makedirs(out)
        with open(os.path.join(out, "search_provenance.json"), "w") as fh:
            json.dump(prov, fh, indent=2)
        return out

    def _record(self, out):
        rec = mm.search_record(None, os.path.join(out, "search_provenance.json"), None)
        rec["mbr"] = {"value": True, "label": "on", "source": "dda.cfg"}   # --reanalyse in cfg
        return rec

    def test_the_run_registry_says_first_pass_not_match_between_runs(self):
        import record_run
        with tempfile.TemporaryDirectory() as d:
            out = self._out(d)
            rec = self._record(out)
            text = "\n".join(record_run.render_parameters({"parameters": {"record": rec}}))
            self.assertNotIn("| Match-between-runs |", text)
            self.assertIn("| Report quantified |", text)
            self.assertIn(dda_first_pass.CAVEAT, text)
            s = {"search": {"out_dir": out, "parameters": {"record": rec}, "status": "completed"},
                 "prot": {"prot": "PROT_9001"}}
            hits = [n for n in record_run.data_quality_notes(s)
                    if n["what"] == dda_first_pass.CAVEAT]
            self.assertEqual(len(hits), 1)
            self.assertEqual(hits[0]["severity"], "WARNING")
            # AUDIT.md said it already: the audit's finding is the one note
            s["analysis"] = {"audit": [{"check": "dda_first_pass", "status": "WARN",
                                        "message": dda_first_pass.CAVEAT}]}
            hits = [n for n in record_run.data_quality_notes(s)
                    if dda_first_pass.CAVEAT in n["what"]]
            self.assertEqual([n["source"] for n in hits], ["AUDIT.md"])
            # an ordinary search keeps its match-between-runs row
            rec2 = dict(rec, dda_first_pass=None)
            text = "\n".join(record_run.render_parameters({"parameters": {"record": rec2}}))
            self.assertIn("| Match-between-runs |", text)

    def test_the_slack_post_and_checkpoint_say_it(self):
        import notify_slack
        with tempfile.TemporaryDirectory() as d:
            out = self._out(d)
            self.assertEqual(dda_first_pass.caution_for(out), dda_first_pass.CAVEAT)
            f = notify_slack.search_facts(out, exit_code=0)
            self.assertEqual(f["dda_first_pass"], dda_first_pass.CAVEAT)
            msg = notify_slack.render(f)
            self.assertIn("first-pass report", msg["text"])
            self.assertNotIn("fell back, not measured", msg["text"])
            self.assertTrue(any(dda_first_pass.CAVEAT.replace("'", "") in
                                json.dumps(b).replace("'", "") for b in msg["blocks"]))
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
            self.assertIn(dda_first_pass.CAVEAT, st.get("warnings") or [], st)
        with tempfile.TemporaryDirectory() as d:          # an ordinary search: nothing said
            out = self._out(d, {"engine": "diann"})
            self.assertIsNone(dda_first_pass.caution_for(out))
            self.assertIsNone(notify_slack.search_facts(out, exit_code=0).get("dda_first_pass"))


class Watcher(unittest.TestCase):
    def watch(self, prov, mode="slurm"):
        with tempfile.TemporaryDirectory() as d, tempfile.TemporaryDirectory() as o:
            if prov is not None:
                with open(os.path.join(o, "search_provenance.json"), "w") as fh:
                    json.dump(prov, fh, indent=2)
            sacct = 'echo "  COMPLETED"'
            for name, body in (("sacct", sacct), ("squeue", "")):
                p = os.path.join(d, name)
                with open(p, "w") as fh:
                    fh.write(f"#!/bin/bash\n{body}\n")
                os.chmod(p, 0o755)
            if mode == "chain":
                with open(os.path.join(o, "jobs.txt"), "w") as fh:
                    fh.write("201\n202\n")
                args = ["bash", WATCH, "--all", o]
            else:
                args = ["bash", WATCH, "--slurm", "999", "--out", o]
            env = dict(os.environ, PATH=d + os.pathsep + os.environ.get("PATH", ""))
            r = subprocess.run(args, capture_output=True, text=True, env=env, timeout=60)
            self.assertEqual(r.returncode, 0, r.stderr)
            return json.loads(r.stdout)

    def test_watcher_reports_the_fallback_as_done_with_caution(self):
        for mode in ("slurm", "chain"):
            with self.subTest(mode=mode):
                out = self.watch(FALLBACK_PROV, mode)
                self.assertTrue(out["done"])
                self.assertFalse(out["failed"])
                self.assertIs(out["dda_first_pass"], True)
                self.assertEqual(out["caution"], dda_first_pass.caution(FALLBACK_PROV))

    def test_watcher_old_provenance_has_no_caution(self):
        """Review Focus 1: a provenance from before 2.14.2 has no report_used."""
        old = {k: v for k, v in FALLBACK_PROV.items() if k not in ("report_used", "dda_first_pass")}
        for mode in ("slurm", "chain"):
            for prov in (old, None):
                with self.subTest(mode=mode, prov=bool(prov)):
                    out = self.watch(prov, mode)
                    self.assertTrue(out["done"])
                    self.assertNotIn("dda_first_pass", out)
                    self.assertNotIn("caution", out)


if __name__ == "__main__":
    unittest.main()
