#!/usr/bin/env python3
"""
The 5-step chain's two passes, for every chain -- DIA as well as DDA (the proteomics review of
fix/2.9-dda-chain, 2026-09-29):

  F1  the final pass can LOSE what the first pass found, and nothing compared them. SET28 (6
      Exploris DDA hair runs, DIA-NN 2.7.0), precursors per run, step 3 -> step 5: 103 -> 10,
      566 -> 529, 126 -> 0, 138 -> 7, 269 -> 183, 165 -> 1. Step 5 now tabulates both passes
      (pass_comparison.py), warns on runs that fell by more than half, records the table in
      search_provenance.json, and audit_results.py puts it in AUDIT.md. The Methods describe the
      first-pass report when the DE used it.
  F2  the array tasks had no --out, so every task wrote report.parquet & co. into <out> itself.
  F8  the Methods said "with match-between-runs" from --reanalyse in the cfg -- which the chain
      strips -- instead of describing the chain's own two passes from search_mode.

Fixture stats files carry the SET28 numbers under generic run names. Hermetic: nothing here
reads a real run, a real provenance or the network. The rows-after-q-filter checks build small
parquet reports and run only where pyarrow is installed (CI installs none).
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:                       # the stdlib-only CI: those tests are skipped
    pa = pq = None

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

from job_env import job_env  # noqa: E402  (env for anything that runs generated bash)

import diann_parallel as dp  # noqa: E402
import make_methods as mm  # noqa: E402
import pass_comparison as pc  # noqa: E402
import record_run  # noqa: E402

COMPARE = os.path.join(SCRIPTS, "pass_comparison.py")
AUDIT = os.path.join(SCRIPTS, "audit_results.py")
GENERATE = os.path.join(SCRIPTS, "diann_parallel.py")

# SET28, step 3 -> step 5: (run, first precursors, first PGs, final precursors, final PGs)
SET28 = [("run1", 103, 0, 10, 0), ("run2", 566, 184, 529, 0), ("run3", 126, 0, 0, 0),
        ("run4", 138, 0, 7, 0), ("run5", 269, 116, 183, 106), ("run6", 165, 0, 1, 0)]
FLAGGED = ["run1", "run3", "run4", "run6"]


def read(path):
    with open(path) as fh:
        return fh.read()


def write_stats(path, rows):
    """A DIA-NN <report>.stats.tsv: File.Name (a raw path), then the per-run counts."""
    with open(path, "w") as fh:
        fh.write("File.Name\tPrecursors.Identified\tProteins.Identified\tTotal.Quantity\n")
        for run, pr, pg in rows:
            fh.write(f"/data/raw/{run}.raw\t{pr}\t{pg}\t1e+08\n")


def chain_out(d, rows=SET28):
    """<out> with both passes' stats files, as steps 3 and 5 leave them."""
    out = os.path.join(d, "out")
    os.makedirs(out)
    write_stats(os.path.join(out, "step3_assembly.stats.tsv"), [(r, a, b) for r, a, b, _, _ in rows])
    write_stats(os.path.join(out, "report.stats.tsv"), [(r, c, e) for r, _, _, c, e in rows])
    return out, os.path.join(out, "step3_assembly.parquet"), os.path.join(out, "report.parquet")


def cohort(d, n=6):
    paths = []
    for i in range(n):
        p = os.path.join(d, f"run{i}.mzML")
        with open(p, "w") as fh:
            fh.write("<mzML/>")
        paths.append(p)
    return paths


def generate(d, cfg_text, out_name="gen"):
    cfg = os.path.join(d, f"{out_name}.cfg")
    with open(cfg, "w") as fh:
        fh.write(cfg_text)
    fasta = os.path.join(d, "db.fasta")
    with open(fasta, "w") as fh:
        fh.write(">sp|P1|X\nPEPTIDER\n")
    out = os.path.join(d, out_name)
    p = subprocess.run([sys.executable, GENERATE, "--diann", "/bin/true", "--raw", *cohort(d),
                        "--fasta", fasta, "--out", out, "--cfg", cfg],
                       capture_output=True, text=True, timeout=120)
    return p, out


PINNED = "--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n--window 7\n"

# Rows passing run_de.R's q filter, per run, SET28-shaped: the four flagged runs have none in
# EITHER report; only run5 -- not flagged -- has more in the first pass (SET28: 396 vs 147).
Q_ROWS = {"first": {"run1": 0, "run2": 30, "run3": 0, "run4": 0, "run5": 40, "run6": 0},
          "final": {"run1": 0, "run2": 30, "run3": 0, "run4": 0, "run5": 15, "run6": 0}}
Q_COLS = ["Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "PG.Q.Value", "Global.Q.Value",
          "Global.PG.Q.Value"]


def write_report(path, passing):
    """A report.parquet: for each run, `passing` rows that pass every q-column, plus 5 rows that
    fail one (Q.Value 0.02 > 0.01) -- so only the filter separates them."""
    runs, cols = [], {c: [] for c in Q_COLS}
    for run, n in passing.items():
        for i in range(n + 5):
            runs.append(run)
            for c in Q_COLS:
                cols[c].append(0.02 if (c == "Q.Value" and i >= n) else 0.001)
    pq.write_table(pa.table({"Run": runs, **cols}), path)


def with_reports(out, first=Q_ROWS["first"], final=Q_ROWS["final"]):
    write_report(os.path.join(out, "step3_assembly.parquet"), first)
    write_report(os.path.join(out, "report.parquet"), final)


# ------------------------------------------------------------------------------------------
# F1: the comparison itself
# ------------------------------------------------------------------------------------------
class PassComparison(unittest.TestCase):
    def test_the_set28_numbers_flag_the_runs_that_lost_more_than_half(self):
        with tempfile.TemporaryDirectory() as d:
            _, first, final = chain_out(d)
            rec = pc.compare(first, final)
        self.assertEqual(rec["flagged"], FLAGGED)
        self.assertEqual((rec["n_runs"], rec["n_flagged"]), (6, 4))
        by = {r["run"]: r for r in rec["runs"]}
        self.assertEqual(by["run3"]["precursor_change"], -1.0)
        self.assertEqual(by["run2"]["first"],
                         {"precursors": 566, "proteins": 184, "rows_after_q": None})
        self.assertEqual(by["run2"]["final"],
                         {"precursors": 529, "proteins": 0, "rows_after_q": None})
        self.assertFalse(by["run5"]["flagged"], "269 -> 183 kept more than half")
        self.assertIn("step3_assembly.parquet", rec["advice"])
        # no parquet here: the post-filter rows are not counted, and nothing is recommended
        self.assertFalse(rec["post_filter"]["available"])
        self.assertIsNone(rec["switch_recommended_for"])
        self.assertIn("before protein-group FDR", rec["advice"])
        self.assertIn("could not be counted", rec["advice"])
        self.assertNotIn("offer the first-pass", rec["advice"])
        self.assertEqual(pc.table(rec).count("**FLAGGED**"), 4)
        self.assertEqual(len(pc.warnings(rec)), 4)

    def test_a_final_pass_that_holds_up_flags_nothing(self):
        rows = [("a", 5000, 800, 6100, 850), ("b", 4800, 790, 2500, 700)]
        with tempfile.TemporaryDirectory() as d:
            _, first, final = chain_out(d, rows)
            rec = pc.compare(first, final)
        self.assertEqual(rec["flagged"], [])
        self.assertEqual(pc.warnings(rec), [])
        self.assertIn("no run lost", rec["advice"])

    def test_a_run_the_first_pass_found_nothing_in_is_not_flagged(self):
        rows = [("a", 0, 0, 0, 0), ("b", 0, 0, 12, 1)]
        with tempfile.TemporaryDirectory() as d:
            _, first, final = chain_out(d, rows)
            rec = pc.compare(first, final)
        self.assertEqual(rec["flagged"], [])
        self.assertIsNone(rec["runs"][0]["precursor_change"])

    def test_the_cli_prints_warns_writes_and_records_in_the_provenance(self):
        with tempfile.TemporaryDirectory() as d:
            out, first, final = chain_out(d)
            prov = os.path.join(out, "search_provenance.json")
            with open(prov, "w") as fh:
                json.dump({"engine": "diann", "search_mode": "parallel_5step"}, fh)
            res = subprocess.run([sys.executable, COMPARE, "--first-report", first,
                                  "--final-report", final, "--out",
                                  os.path.join(out, "pass_comparison.json"), "--provenance", prov],
                                 capture_output=True, text=True, timeout=60)
            self.assertEqual(res.returncode, 0, res.stderr)
            self.assertIn("| run3 | 126 | 0 | -100% **FLAGGED** |", res.stdout)
            self.assertEqual(res.stderr.count("WARNING: run"), 4)
            self.assertIn("before protein-group FDR", res.stderr)
            saved = json.loads(read(os.path.join(out, "pass_comparison.json")))
            recorded = json.loads(read(prov))
            self.assertEqual(recorded["engine"], "diann", "the provenance kept what it had")
            self.assertEqual(recorded["pass_comparison"], saved)
            self.assertFalse(os.path.exists(prov + ".tmp"))

    @unittest.skipUnless(pq, "pyarrow not installed")
    def test_the_rows_after_run_des_filter_decide_the_switch(self):
        """The SET28 shape: four runs flagged, none of them with a single row after the q filter
        in either report; run5, unflagged, has more in the first pass. Only run5 is offered."""
        with tempfile.TemporaryDirectory() as d:
            out, first, final = chain_out(d)
            with_reports(out)
            rec = pc.compare(first, final)
            text = pc.table(rec)
        self.assertTrue(rec["post_filter"]["available"])
        self.assertEqual(rec["flagged"], FLAGGED)
        self.assertEqual(rec["switch_recommended_for"], ["run5"])
        by = {r["run"]: r for r in rec["runs"]}
        self.assertEqual((by["run3"]["first"]["rows_after_q"], by["run3"]["final"]["rows_after_q"]),
                         (0, 0))
        self.assertEqual((by["run5"]["first"]["rows_after_q"], by["run5"]["final"]["rows_after_q"]),
                         (40, 15))
        self.assertIn("more rows for 1 run(s): run5", rec["advice"])
        self.assertIn("| run5 | 269 | 183 | -32% | 116 | 106 | 40 **first pass has more** | 15 |",
                      text)
        self.assertIn("before protein-group FDR", text)
        self.assertEqual(rec["post_filter"]["cutoffs"]["PG.Q.Value"], 0.05)
        self.assertEqual(rec["post_filter"]["cutoffs"]["Q.Value"], 0.01)

    @unittest.skipUnless(pq, "pyarrow not installed")
    def test_flagged_runs_the_first_pass_cannot_recover_are_not_offered(self):
        with tempfile.TemporaryDirectory() as d:
            out, first, final = chain_out(d)
            same = dict(Q_ROWS["final"])
            with_reports(out, first=same, final=same)
            rec = pc.compare(first, final)
        self.assertEqual(rec["switch_recommended_for"], [])
        self.assertIn("would not recover them", rec["advice"])
        self.assertIn("keep the final report", rec["advice"])

    @unittest.skipUnless(pq, "pyarrow not installed")
    def test_the_filter_is_run_des_as_limpa_applies_it(self):
        """A row goes only when a value EXCEEDS its cutoff (limpa: which(x > cutoff)), a missing
        value stays, PG.Q.Value is held at 0.05 (cutoff_for), and a column the report lacks is
        not filtered on (fdr_columns)."""
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "r.parquet")
            pq.write_table(pa.table({
                "Run": ["a", "a", "a", "a", "a"],
                "Q.Value": [0.01, 0.0100001, None, 0.001, 0.001],
                "Lib.Q.Value": [0.0, 0.0, 0.0, 0.0, 0.0],
                "Lib.PG.Q.Value": [0.0, 0.0, 0.0, 0.0, 0.0],
                "PG.Q.Value": [0.0, 0.0, 0.0, 0.03, 0.06]}), path)
            rows, cols, why = pc.rows_after_q(path)
        self.assertIsNone(why)
        self.assertEqual(cols, ["Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "PG.Q.Value"])
        # kept: == cutoff, missing, PG 0.03 <= 0.05; dropped: 0.0100001, PG 0.06
        self.assertEqual(rows, {"a": 3})

    def test_a_missing_stats_file_is_a_stated_failure(self):
        with tempfile.TemporaryDirectory() as d:
            out, first, final = chain_out(d)
            os.remove(os.path.join(out, "step3_assembly.stats.tsv"))
            res = subprocess.run([sys.executable, COMPARE, "--first-report", first,
                                  "--final-report", final], capture_output=True, text=True,
                                 timeout=60)
        self.assertEqual(res.returncode, 1)
        self.assertIn("step3_assembly.stats.tsv", res.stderr)


# ------------------------------------------------------------------------------------------
# F1 + F2: what the chain generates
# ------------------------------------------------------------------------------------------
class ChainScripts(unittest.TestCase):
    def test_step5_compares_the_passes_and_step3_clears_a_stale_first_pass(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = generate(d, PINNED)
            self.assertEqual(p.returncode, 0, p.stderr)
            s3 = read(os.path.join(out, "step3_assembly.sbatch"))
            s5 = read(os.path.join(out, "step5_report.sbatch"))
            info = json.loads(p.stdout)
        first = os.path.join(out, dp.FIRST_PASS_REPORT)
        self.assertIn(f'"{first}" "{os.path.join(out, "step3_assembly.stats.tsv")}"', s3)
        self.assertIn(f'--out "{first}"', s3)   # quoted: a folder may hold a space
        line = next(ln for ln in s5.splitlines() if "pass_comparison.py" in ln)
        self.assertIn(f"--first-report {first}", line)
        self.assertIn(f"--final-report {os.path.join(out, 'report.parquet')}", line)
        self.assertIn(f"--provenance {os.path.join(out, 'search_provenance.json')}", line)
        self.assertIn('|| echo "WARNING', line, "a comparison that cannot be made is said")
        self.assertLess(s5.index("OK: report built"), s5.index("pass_comparison.py"))
        self.assertEqual(info["pass_comparison"]["produced"], "runtime")
        self.assertEqual(info["pass_comparison"]["first_pass_report"], first)

    def test_the_step5_command_runs_as_generated(self):
        """The generated line itself, in bash, against fixture stats: quoting and paths."""
        with tempfile.TemporaryDirectory() as d:
            p, out = generate(d, PINNED)
            self.assertEqual(p.returncode, 0, p.stderr)
            line = next(ln for ln in read(os.path.join(out, "step5_report.sbatch")).splitlines()
                        if "pass_comparison.py" in ln)
            for name, col in (("step3_assembly.stats.tsv", 1), ("report.stats.tsv", 3)):
                write_stats(os.path.join(out, name), [(r[0], r[col], r[col + 1]) for r in SET28])
            with open(os.path.join(out, "search_provenance.json"), "w") as fh:
                json.dump({"search_mode": "parallel_5step"}, fh)
            res = subprocess.run(["bash", "-c", line], capture_output=True, text=True,
                                 env=job_env(d), timeout=60)
            self.assertEqual(res.returncode, 0, res.stderr)
            self.assertIn("**FLAGGED**", res.stdout)
            self.assertEqual(json.loads(read(os.path.join(out, "search_provenance.json")))
                             ["pass_comparison"]["flagged"], FLAGGED)

    def test_array_tasks_write_their_reports_into_their_own_folders(self):
        """F2: without --out every task wrote report.parquet, report.stats.tsv,
        report-lib.parquet and report.log.txt into <out>, concurrently."""
        with tempfile.TemporaryDirectory() as d:
            p, out = generate(d, PINNED, "noxic")
            self.assertEqual(p.returncode, 0, p.stderr)
            s2 = read(os.path.join(out, "step2_firstpass.sbatch"))
            s4 = read(os.path.join(out, "step4_finalpass.sbatch"))
            p, xout = generate(d, PINNED + "--xic 10\n", "xic")
            self.assertEqual(p.returncode, 0, p.stderr)
            s4x = read(os.path.join(xout, "step4_finalpass.sbatch"))
            s3 = read(os.path.join(out, "step3_assembly.sbatch"))
        task = "t${SLURM_ARRAY_TASK_ID}.parquet"
        self.assertIn(f'mkdir -p "{out}/firstpass"', s2)
        self.assertIn(f'--out "{out}/firstpass/{task}"', s2)
        self.assertIn(f'mkdir -p "{out}/finalpass"', s4)
        self.assertIn(f'--out "{out}/finalpass/{task}"', s4)
        self.assertIn(f'--out "{xout}/xic/{task}"', s4x)
        self.assertNotIn(f"{xout}/finalpass", s4x)
        # steps 3 and 5 still read the per-file .quant folders, not the task reports
        self.assertIn("--use-quant", s3)
        self.assertIn(f'--temp "{out}/quant_step2"', s3)


# ------------------------------------------------------------------------------------------
# F1: AUDIT.md
# ------------------------------------------------------------------------------------------
class AuditReportsThePasses(unittest.TestCase):
    def audit(self, d, de_input, mode, comparison=True, reports=False):
        """AUDIT.md for a DE whose input is <out>/<de_input>, beside a search_provenance.json
        of `mode` -- with step 5's pass_comparison in it, or without."""
        out, first, final = chain_out(d)
        if reports:
            with_reports(out)
        prov = {"search_mode": mode}
        if comparison:
            prov["pass_comparison"] = pc.compare(first, final)
        with open(os.path.join(out, "search_provenance.json"), "w") as fh:
            json.dump(prov, fh)
        tables = os.path.join(d, "tables")
        os.makedirs(tables)
        with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
            json.dump({"input": os.path.join(out, de_input), "adjp": 0.05, "logfc": 1}, fh)
        res = subprocess.run([sys.executable, AUDIT, "--out", "AUDIT.md", "--de-dir", tables],
                             capture_output=True, text=True, timeout=120, cwd=d)
        self.assertEqual(res.returncode, 0, res.stderr)
        found = [f for f in json.loads(read(os.path.join(d, "AUDIT.json")))["findings"]
                 if f["check"] == "first_vs_final_pass"]
        return found, read(os.path.join(d, "AUDIT.md"))

    def test_flagged_runs_are_a_warning_in_audit_md(self):
        with tempfile.TemporaryDirectory() as d:
            found, md = self.audit(d, "report.parquet", "parallel_5step")
        self.assertEqual([f["status"] for f in found], ["WARN"])
        self.assertIn("4 of 6 runs", found[0]["message"])
        self.assertIn("run3 126 -> 0", found[0]["message"])
        self.assertIn("counts before protein-group FDR; compare the rows that pass the q-value "
                      "filter before switching", found[0]["message"])
        self.assertIn("were not counted", found[0]["message"])
        self.assertIn("first_vs_final_pass", md)
        self.assertIn(dp.FIRST_PASS_REPORT, md)

    @unittest.skipUnless(pq, "pyarrow not installed")
    def test_audit_offers_the_first_pass_only_where_it_has_more_rows(self):
        with tempfile.TemporaryDirectory() as d:
            found, _ = self.audit(d, "report.parquet", "parallel_5step", reports=True)
        self.assertEqual([f["status"] for f in found], ["WARN"])
        msg = found[0]["message"]
        self.assertIn("run3 126 -> 0 (0 vs 0 rows after the q filter)", msg)
        self.assertIn("the first pass has more rows for run5", msg)

    def test_a_de_on_the_first_pass_is_noted_not_warned(self):
        with tempfile.TemporaryDirectory() as d:
            found, _ = self.audit(d, dp.FIRST_PASS_REPORT, "parallel_5step")
        self.assertEqual([f["status"] for f in found], ["INFO"])
        self.assertIn("The DE used the first-pass report", found[0]["message"])

    def test_a_chain_with_no_comparison_says_so_and_single_shot_says_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            found, _ = self.audit(d, "report.parquet", "parallel_5step", comparison=False)
        self.assertEqual([f["status"] for f in found], ["INFO"])
        self.assertIn("not in search_provenance.json", found[0]["message"])
        with tempfile.TemporaryDirectory() as d:
            found, _ = self.audit(d, "report.parquet", "single_shot", comparison=False)
        self.assertEqual(found, [])


# ------------------------------------------------------------------------------------------
# F8 + F1: the Methods
# ------------------------------------------------------------------------------------------
CHAIN_WORDS = ("Match-between-runs used DIA-NN's two-pass procedure, run as separate jobs: each "
               "run was first searched against the predicted library; precursors identified "
               "across the experiment were assembled into an experiment-specific empirical "
               "spectral library; every run was then searched again against it. This differs "
               "from the retention-time-alignment-based match-between-runs of DDA software such "
               "as MaxQuant.")


class MethodsDescribeTheSearchThatRan(unittest.TestCase):
    CFG = ("--qvalue 0.01\n--dda\n--fasta-search\n--gen-spec-lib\n--predictor\n--reanalyse\n"
           "--cut K*,R*\n--mass-acc 20\n--mass-acc-ms1 10\n")

    def paragraph(self, d, mode, de_input=None, pass_cmp=None):
        cfg = os.path.join(d, "params.cfg")
        with open(cfg, "w") as fh:
            fh.write(self.CFG)
        prov = {"engine": "diann", "version": "2.7.0", "search_mode": mode, "params_file": cfg,
                "result": {}}
        if pass_cmp:
            prov["pass_comparison"] = pass_cmp
        path = os.path.join(d, "search_provenance.json")
        with open(path, "w") as fh:
            json.dump(prov, fh)
        rec = mm.search_record(cfg, path)
        return rec, mm.search_paragraph(rec, {"input": de_input} if de_input else None)

    def test_the_chain_is_described_as_its_two_passes_not_as_diann_mbr(self):
        with tempfile.TemporaryDirectory() as d:
            rec, text = self.paragraph(d, "parallel_5step")
        self.assertIn(CHAIN_WORDS, text)
        self.assertNotIn("with match-between-runs.", text)
        self.assertNotIn(") (", text, "back-to-back parentheticals")
        self.assertIn("DIA-NN searched the spectra in its DDA mode (--dda), which it describes as "
                      "beta-stage support.", text)
        table = "\n".join(record_run.render_parameters({"parameters": {"record": rec,
                                                                         "why": {}}}))
        self.assertIn("| Match-between-runs | two-pass, as separate jobs", table)

    def test_a_single_shot_search_keeps_diann_mbr(self):
        with tempfile.TemporaryDirectory() as d:
            _, text = self.paragraph(d, "single_shot")
        self.assertIn("predictor), with match-between-runs.", text)
        self.assertNotIn("two-pass procedure", text)

    def test_the_first_pass_deliverable_is_described_as_such(self):
        with tempfile.TemporaryDirectory() as d:
            _, text = self.paragraph(d, "parallel_5step",
                                     de_input=os.path.join(d, dp.FIRST_PASS_REPORT),
                                     pass_cmp={"n_runs": 6, "n_flagged": 4, "runs": []})
        self.assertIn("Quantification used the experiment-wide first-pass report: each run was "
                      "searched once against the predicted library, with no match-between-runs",
                      text)
        # the fact, not a recovery claim: the flagged runs may have nothing after the q filter
        self.assertNotIn("because", text)
        self.assertNotIn("retained less than half", text)
        self.assertNotIn("two-pass procedure", text)
        self.assertNotIn("with match-between-runs.", text)


if __name__ == "__main__":
    unittest.main()
