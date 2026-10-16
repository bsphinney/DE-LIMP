#!/usr/bin/env python3
"""collect_conditions.py --validate --exclude: a run left out of the DE on purpose is a recorded
decision, not a failed check.

PROT_0000 (2026-10-08; skill issue record 2026-10-08_<user>_PROT_0000+….md): the search had
13 runs, two of them re-runs, and the DE analysed 11. conditions.csv listed the 11, so `--validate
--against report.parquet` exited 1 ("2 report runs missing from metadata"). run_de.R already
analyses only the metadata's runs, and step 8 says the validation must pass, so the agent went
past a failed check on its own judgment and the reason was written only in decisions.md.

Guards:
  * --exclude <run> --reason "<why>" (repeatable): the PROT_0000 case validates, and the runs and
    their reasons are written to <csv>.decisions.json (excluded_runs), keeping the record's other
    answers;
  * a run missing WITHOUT --exclude still fails, with the same problem as before;
  * every --exclude needs its own --reason; an excluded run must be in the report and must have no
    metadata row; --exclude needs --against;
  * run_de.R copies the reasons into de_provenance.json (runs_left_out) and methods.txt, and the
    Methods (make_methods.de_runs_left_out_sentence) and the report's callout
    (make_analysis_html.runs_left_out_note) list each run with its reason -- "not recorded" for a
    run left out with none.
NOTE: CI runs Python only and skips the R-invoking class here -- run it locally.
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import collect_conditions as cc  # noqa: E402
import make_analysis_html as mah  # noqa: E402
import make_methods as mm  # noqa: E402
from test_run_de_block import r_has, rscript  # noqa: E402

SCRIPT = os.path.join(SCRIPTS, "collect_conditions.py")
RUN_DE = os.path.join(SCRIPTS, "run_de.R")

# PROT_0000's shape: 11 samples analysed, two of them injected again (13 runs searched).
ANALYSED = [f"RUN_{i:02d}" for i in range(1, 12)]
RERUNS = ["RUN_03_rerun", "RUN_07_rerun"]
REASONS = {"RUN_03_rerun": "re-run of RUN_03 after a pressure fault; the first injection is analysed",
           "RUN_07_rerun": "re-run of RUN_07; the first injection is analysed"}


def write_report(tmp, runs):
    """A DIA-NN-style TSV report: collect_conditions reads only its Run column."""
    path = os.path.join(tmp, "report.tsv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Run", "Precursor.Id"])
        for r in runs:
            w.writerow([r, "PEPTIDEK2"])
    return path


def write_conditions(tmp, runs, name="conditions.csv"):
    path = os.path.join(tmp, name)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File.Name", "Group"])
        for i, r in enumerate(runs):
            w.writerow([r, "Control" if i < len(runs) // 2 else "Treated"])
    return path


def validate(csv_path, report=None, *extra):
    args = [sys.executable, SCRIPT, "--validate", csv_path]
    if report:
        args += ["--against", report]
    return subprocess.run(args + list(extra), capture_output=True, text=True)


def exclude_args(reasons):
    out = []
    for run, why in reasons.items():
        out += ["--exclude", run, "--reason", why]
    return out


class ExcludedRunsValidate(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.report = write_report(self.tmp, ANALYSED + RERUNS)
        self.csv = write_conditions(self.tmp, ANALYSED)
        self.dec = self.csv + ".decisions.json"

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_prot_0803_validates_and_records_the_reasons(self):
        p = validate(self.csv, self.report, *exclude_args(REASONS))
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        out = json.loads(p.stdout)
        self.assertTrue(out["valid"])
        self.assertEqual(out["problems"], [])
        self.assertEqual((out["n_report_runs"], out["n_analysed"]), (13, 11))
        self.assertEqual({x["run"]: x["reason"] for x in out["excluded_runs"]}, REASONS)
        self.assertEqual(out["decisions_file"], os.path.abspath(self.dec))
        with open(self.dec, encoding="utf-8") as fh:
            rec = json.load(fh)
        ex = rec[cc.EXCLUDED_KEY]
        self.assertEqual({x["run"]: x["reason"] for x in ex["runs"]}, REASONS)
        self.assertEqual((ex["n_report_runs"], ex["n_analysed"]), (13, 11))
        self.assertEqual(ex["against"], os.path.abspath(self.report))
        self.assertEqual(rec["conditions_csv"], os.path.abspath(self.csv))

    def test_a_missing_run_without_exclude_still_fails_as_before(self):
        p = validate(self.csv, self.report)
        self.assertEqual(p.returncode, 1)
        out = json.loads(p.stdout)
        self.assertFalse(out["valid"])
        self.assertIn(f"2 report runs missing from metadata: {sorted(RERUNS)[:5]}...", out["problems"])
        self.assertIn("--exclude", out["fix"])       # how a deliberate exclusion is recorded
        self.assertFalse(os.path.exists(self.dec))
        # one of the two excluded: the other still fails, by name
        p = validate(self.csv, self.report, "--exclude", RERUNS[0], "--reason", REASONS[RERUNS[0]])
        self.assertEqual(p.returncode, 1)
        out = json.loads(p.stdout)
        self.assertIn(f"1 report runs missing from metadata: {[RERUNS[1]]}...", out["problems"])
        self.assertFalse(os.path.exists(self.dec), "a failed validation recorded an exclusion")

    def test_every_exclude_needs_its_own_reason(self):
        p = validate(self.csv, self.report, "--exclude", RERUNS[0])
        self.assertEqual(p.returncode, 2)
        self.assertIn(f"--exclude {RERUNS[0]} has no --reason", p.stderr)
        # one reason for two runs belongs to the second: the first has none
        p = validate(self.csv, self.report, "--exclude", RERUNS[0], "--exclude", RERUNS[1],
                     "--reason", "re-runs")
        self.assertEqual(p.returncode, 2)
        self.assertIn(f"--exclude {RERUNS[0]} has no --reason", p.stderr)
        p = validate(self.csv, self.report, "--exclude", RERUNS[0], "--reason", "  ")
        self.assertEqual(p.returncode, 2)
        p = validate(self.csv, self.report, "--reason", "why")
        self.assertEqual(p.returncode, 2)
        self.assertIn("--reason must follow the --exclude", p.stderr)
        self.assertFalse(os.path.exists(self.dec))

    def test_an_excluded_run_the_report_does_not_have_errors(self):
        p = validate(self.csv, self.report, *exclude_args(REASONS),
                     "--exclude", "RUN_99_rerun", "--reason", "typo")
        self.assertEqual(p.returncode, 1)
        out = json.loads(p.stdout)
        self.assertFalse(out["valid"])
        self.assertTrue(any("the report does not have: ['RUN_99_rerun']" in x for x in out["problems"]),
                        out["problems"])
        self.assertFalse(os.path.exists(self.dec))

    def test_an_excluded_run_with_a_metadata_row_errors(self):
        both = write_conditions(self.tmp, ANALYSED + RERUNS, name="all.csv")
        p = validate(both, self.report, *exclude_args(REASONS))
        self.assertEqual(p.returncode, 1)
        out = json.loads(p.stdout)
        self.assertTrue(any("still have a metadata row" in x for x in out["problems"]), out["problems"])

    def test_exclude_needs_validate_against_a_report(self):
        p = validate(self.csv, None, *exclude_args(REASONS))
        self.assertEqual(p.returncode, 2)
        self.assertIn("--exclude needs --against", p.stderr)
        p = subprocess.run([sys.executable, SCRIPT, "--list-runs", "--runs", "a,b",
                            "--exclude", "a", "--reason", "x"], capture_output=True, text=True)
        self.assertEqual(p.returncode, 2)
        p = validate(self.csv, self.report, "--exclude", RERUNS[0], "--reason", "a",
                     "--exclude", RERUNS[0], "--reason", "b")
        self.assertEqual(p.returncode, 2)
        self.assertIn("more than once", p.stderr)

    def test_the_record_keeps_the_other_answers_and_drops_a_stale_exclusion(self):
        with open(self.dec, "w", encoding="utf-8") as fh:
            json.dump({"confirmed_multi_match": {"ctrl": ["a", "b"]},
                       "conditions_csv": os.path.abspath(self.csv)}, fh)
        p = validate(self.csv, self.report, *exclude_args(REASONS))
        self.assertEqual(p.returncode, 0, p.stdout)
        with open(self.dec, encoding="utf-8") as fh:
            rec = json.load(fh)
        self.assertEqual(rec["confirmed_multi_match"], {"ctrl": ["a", "b"]})
        self.assertIn(cc.EXCLUDED_KEY, rec)
        # the re-runs are analysed after all: the exclusion no longer describes this design
        write_conditions(self.tmp, ANALYSED + RERUNS)
        p = validate(self.csv, self.report)
        self.assertEqual(p.returncode, 0, p.stdout)
        with open(self.dec, encoding="utf-8") as fh:
            rec = json.load(fh)
        self.assertNotIn(cc.EXCLUDED_KEY, rec)
        self.assertEqual(rec["confirmed_multi_match"], {"ctrl": ["a", "b"]})

    def test_a_record_holding_only_the_exclusion_is_removed_with_it(self):
        self.assertEqual(validate(self.csv, self.report, *exclude_args(REASONS)).returncode, 0)
        self.assertTrue(os.path.exists(self.dec))
        write_conditions(self.tmp, ANALYSED + RERUNS)
        self.assertEqual(validate(self.csv, self.report).returncode, 0)
        self.assertFalse(os.path.exists(self.dec))


def left_out_prov(reasons):
    return {"display_label": "DPC-Quant + limma", "runs_left_out": {
        "determined": True, "n_report_runs": 13, "n_analysed": 11,
        "runs": [{"run": r, "reason": reasons.get(r)} for r in RERUNS]}}


class TheMethodsAndTheReportSaySo(unittest.TestCase):
    def test_the_methods_list_each_run_with_its_reason(self):
        para = mm.de_paragraph(left_out_prov(REASONS))
        self.assertIn("Of the 13 runs searched, 11 were analysed and 2 were left out of the "
                      "differential-expression analysis: RUN_03_rerun (re-run of RUN_03 after a "
                      "pressure fault; the first injection is analysed); RUN_07_rerun (re-run of "
                      "RUN_07; the first injection is analysed).", para)

    def test_a_run_with_no_recorded_reason_says_not_recorded(self):
        s = mm.de_runs_left_out_sentence(left_out_prov({RERUNS[0]: REASONS[RERUNS[0]]}))
        self.assertIn(f"RUN_07_rerun (reason {mm.NOT_RECORDED})", s)

    def test_the_report_callout_uses_the_same_words(self):
        prov = left_out_prov(REASONS)
        n = mah.runs_left_out_note(prov)
        self.assertEqual((n["kind"], n["text"]), ("info", mm.de_runs_left_out_sentence(prov)))
        self.assertEqual(mah.runs_left_out_note(left_out_prov({}))["kind"], "warning")

    def test_a_record_that_could_not_be_read_says_not_read_not_not_recorded(self):
        """2.12 review: run_de.R records `reasons_unreadable` when the decisions record exists but
        cannot be read. The user DID record reasons, so "not recorded" would be wrong."""
        prov = left_out_prov({})
        prov["runs_left_out"]["reasons_unreadable"] = "/x/conditions.csv.decisions.json could not be read (parse error)"
        s = mm.de_runs_left_out_sentence(prov)
        for run in RERUNS:
            self.assertIn(f"{run} (reason not read: the decisions record could not be read [confirm])", s)
        self.assertEqual(mm.REASON_NOT_READ, "reason not read: the decisions record could not be read [confirm]")
        self.assertNotIn(mm.NOT_RECORDED, s)
        self.assertNotIn("/x/", s, "a path in the Methods")
        self.assertIn(s, mm.de_paragraph(prov))
        n = mah.runs_left_out_note(prov)
        self.assertEqual((n["kind"], n["text"]), ("warning", s))
        # a run whose reason WAS read keeps it; the flag only fills the missing ones
        prov = left_out_prov({RERUNS[0]: REASONS[RERUNS[0]]})
        prov["runs_left_out"]["reasons_unreadable"] = "x"
        s = mm.de_runs_left_out_sentence(prov)
        self.assertIn(f"{RERUNS[0]} ({REASONS[RERUNS[0]]})", s)
        self.assertIn(f"{RERUNS[1]} ({mm.REASON_NOT_READ})", s)
        self.assertEqual(mah.runs_left_out_note(prov)["kind"], "warning")

    def test_nothing_is_said_when_every_run_was_analysed_or_the_record_is_older(self):
        none_left = {"runs_left_out": {"determined": True, "n_report_runs": 11, "n_analysed": 11,
                                       "runs": []}}
        for prov in ({}, none_left, {"runs_left_out": {"determined": False, "note": "x"}}):
            self.assertIsNone(mm.de_runs_left_out_sentence(prov))
            self.assertIsNone(mah.runs_left_out_note(prov))
            self.assertNotIn("left out", mm.de_paragraph(dict(prov, display_label="X")))

    def test_run_de_reads_the_key_and_file_collect_conditions_writes(self):
        with open(RUN_DE) as fh:
            src = fh.read()
        self.assertIn("$" + cc.EXCLUDED_KEY, src)
        self.assertTrue(cc.decisions_path("x.csv").startswith("x.csv"))
        self.assertIn(f'paste0(meta_path, "{cc.decisions_path("")}")', src)


# 8 runs: A1-A3 and B1-B3 analysed, A2 and B2 injected again (the two re-runs left out). Low
# precursors go undetected (limpa's detection-probability curve needs missing values to fit).
SYNTH_R = r'''
set.seed(3)
runs <- c("A1", "A2", "A3", "B1", "B2", "B3", "A2_rerun", "B2_rerun")
grp <- c("A", "A", "A", "B", "B", "B", "A", "B")
nr <- length(runs)
out <- list(); pgq <- list()
for (i in 1:150) {
  pg <- sprintf("P%05d", i); base <- runif(1, 12, 20)
  eff <- if (i <= 20) 0.8 else if (i <= 40) -0.8 else 0
  run_lv <- base + (grp == "B") * eff + rnorm(nr, 0, 0.25)
  pgq[[i]] <- data.frame(Protein.Group = pg, Run = runs, PG.MaxLFQ = 2^run_lv)
  for (j in 1:4) {
    lv <- run_lv + rnorm(1, 0, 1)
    int <- 2^(lv + rnorm(nr, 0, 0.1))
    keep <- runif(nr) >= plogis(-(lv - 12) * 2)
    if (!any(keep)) next
    out[[length(out) + 1]] <- data.frame(Run = runs[keep], Precursor.Id = sprintf("PEP%dK%d2", i, j),
      Protein.Group = pg, Protein.Ids = pg, Protein.Names = paste0("G", i, "_X"),
      Genes = paste0("G", i), Proteotypic = 1L, Precursor.Normalised = int[keep],
      Precursor.Quantity = int[keep], Q.Value = 0.001, Lib.Q.Value = 0.001,
      Lib.PG.Q.Value = 0.001, PG.Q.Value = 0.001, Global.Q.Value = 0.001,
      Global.PG.Q.Value = 0.001, stringsAsFactors = FALSE)
  }
}
d <- merge(do.call(rbind, out), do.call(rbind, pgq), by = c("Protein.Group", "Run"))
arrow::write_parquet(d, file.path(OUT, "report.parquet"))
'''
E2E_REASONS = {"A2_rerun": "re-run of A2; the first injection is analysed",
               "B2_rerun": "re-run of B2; the first injection is analysed"}


@unittest.skipUnless(r_has("limpa", "limma", "arrow", "dplyr", "tidyr", "jsonlite"),
                     "needs R with limpa/limma/arrow/dplyr/tidyr/jsonlite")
class RunDeRecordsTheExclusion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        analysed = ["A1", "A2", "A3", "B1", "B2", "B3"]
        report_tsv = write_report(cls.tmp, analysed + list(E2E_REASONS))   # the same Run names
        meta = write_conditions(cls.tmp, analysed)
        p = validate(meta, report_tsv, *exclude_args(E2E_REASONS))
        assert p.returncode == 0, p.stdout + p.stderr
        # the same design with no decisions record beside it: reasons are not recorded
        bare = os.path.join(cls.tmp, "bare")
        os.makedirs(bare)
        bare_meta = shutil.copy(meta, os.path.join(bare, "conditions.csv"))
        cls.runs = {}
        for key, m, method in (("dpc", meta, "dpc"), ("ml", meta, "maxlfq"),
                               ("bare", bare_meta, "dpc")):
            out = os.path.join(cls.tmp, "out_" + key)
            p = subprocess.run(["Rscript", RUN_DE, "--input", os.path.join(cls.tmp, "report.parquet"),
                                "--metadata", m, "--outdir", out, "--method", method],
                               capture_output=True, text=True, cwd=cls.tmp)
            cls.runs[key] = (p, out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def prov(self, key):
        p, out = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            return json.load(fh)

    def test_de_provenance_carries_the_reasons_on_both_methods(self):
        for key in ("dpc", "ml"):
            r = self.prov(key)["runs_left_out"]
            self.assertTrue(r["determined"], key)
            self.assertEqual((r["n_report_runs"], r["n_analysed"]), (8, 6), key)
            self.assertEqual({x["run"]: x["reason"] for x in r["runs"]}, E2E_REASONS, key)
            self.assertEqual(self.prov(key)["n_samples"], 6, key)

    def test_methods_txt_and_the_methods_list_them(self):
        _p, out = self.runs["dpc"]
        with open(os.path.join(out, "methods.txt")) as fh:
            text = fh.read()
        self.assertIn("Runs          : 6 of 8 searched runs analysed; left out of the DE:", text)
        self.assertIn("A2_rerun -- re-run of A2; the first injection is analysed", text)
        para = mm.de_paragraph(self.prov("dpc"))
        self.assertIn("Of the 8 runs searched, 6 were analysed and 2 were left out", para)
        self.assertIn("B2_rerun (re-run of B2; the first injection is analysed)", para)

    def test_without_a_record_the_reasons_are_not_recorded_never_invented(self):
        r = self.prov("bare")["runs_left_out"]
        self.assertEqual([x["reason"] for x in r["runs"]], [None, None])
        self.assertIn("--exclude was not used", r["reasons_from"])
        _p, out = self.runs["bare"]
        with open(os.path.join(out, "methods.txt")) as fh:
            self.assertIn("A2_rerun -- reason NOT RECORDED [confirm]", fh.read())

    def test_nothing_left_out_is_an_empty_list(self):
        # every report run analysed: the record says so, and the Methods say nothing
        meta = write_conditions(self.tmp, ["A1", "A2", "A3", "B1", "B2", "B3", "A2_rerun",
                                           "B2_rerun"], name="all.csv")
        out = os.path.join(self.tmp, "out_all")
        p = subprocess.run(["Rscript", RUN_DE, "--input", os.path.join(self.tmp, "report.parquet"),
                            "--metadata", meta, "--outdir", out, "--method", "maxlfq"],
                           capture_output=True, text=True, cwd=self.tmp)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov["runs_left_out"]["runs"], [])
        self.assertIsNone(mm.de_runs_left_out_sentence(prov))


class DecisionsRecordReview212(unittest.TestCase):
    """2.12 review: run_de.R read an unreadable decisions.json as "no reasons" (tryCatch(error =
    NULL)), so every run left out came out NOT RECORDED with nothing said; and collect_conditions
    wrote the record in place, so an interrupted write left exactly such a file."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.meta = write_conditions(self.tmp, ANALYSED)
        self.dec = cc.decisions_path(self.meta)

    def tearDown(self):
        self._tmp.cleanup()

    def leftovers(self):
        return [f for f in os.listdir(self.tmp) if f.endswith(".tmp")]

    # ---- collect_conditions: whole or not at all ---------------------------------------------
    def test_an_interrupted_write_keeps_the_old_record(self):
        from unittest import mock
        old = {"technical_replicates": {"column": "Sample"}, "conditions_csv": self.meta}
        cc.write_decisions(self.dec, old)

        def half(obj, fh, **kw):
            fh.write('{"excluded_runs": {"runs": [')
            raise KeyboardInterrupt
        with mock.patch.object(cc.json, "dump", side_effect=half):
            with self.assertRaises(KeyboardInterrupt):
                cc.write_decisions(self.dec, {"excluded_runs": {}})
        with open(self.dec, encoding="utf-8") as fh:
            self.assertEqual(json.load(fh), old)
        self.assertEqual(self.leftovers(), [])

    def test_validate_exclude_writes_through_the_atomic_writer(self):
        report = write_report(self.tmp, ANALYSED + RERUNS)
        p = validate(self.meta, report, *exclude_args(REASONS))
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        with open(self.dec, encoding="utf-8") as fh:
            rec = json.load(fh)
        self.assertEqual({x["run"]: x["reason"] for x in rec[cc.EXCLUDED_KEY]["runs"]}, REASONS)
        self.assertEqual(self.leftovers(), [])
        with open(SCRIPT, encoding="utf-8") as fh:
            src = fh.read()
        self.assertNotIn('open(dpath, "w"', src, "a decisions record written in place")

    def test_a_link_at_the_temporary_name_is_not_written_through(self):
        theirs = os.path.join(self.tmp, "theirs.txt")
        open(theirs, "w").close()
        os.symlink(theirs, f"{self.dec}.{os.getpid()}.tmp")
        with self.assertRaises(FileExistsError):
            cc.write_decisions(self.dec, {"x": 1})
        with open(theirs) as fh:
            self.assertEqual(fh.read(), "")
        self.assertFalse(os.path.exists(self.dec))

    # ---- run_de.R: an unreadable record is said, on stderr and in the record -----------------
    def record(self, dec_text):
        """runs_left_out_record() from run_de.R (its function definitions only), on two report
        runs the metadata leaves out; stderr is returned too."""
        if dec_text is not None:
            with open(self.dec, "w", encoding="utf-8") as fh:
                fh.write(dec_text)
        out = os.path.join(self.tmp, "rec.json")
        code = f"""
          for (e in parse({json.dumps(RUN_DE)})) if (is.call(e) && identical(e[[1]], as.name('<-')) &&
              is.call(e[[3]]) && identical(e[[3]][[1]], as.name('function')) &&
              as.character(e[[2]]) %in% c('runs_left_out_record', 'runs_left_out_methods_lines')) eval(e)
          r <- runs_left_out_record(c({', '.join(json.dumps(x) for x in ANALYSED + RERUNS)}),
                                    c({', '.join(json.dumps(x) for x in ANALYSED)}), {json.dumps(self.meta)})
          r$methods <- runs_left_out_methods_lines(r)
          jsonlite::write_json(r, {json.dumps(out)}, auto_unbox = TRUE, null = "null", na = "null")"""
        p = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(out) as fh:
            return json.load(fh), p.stderr

    @unittest.skipUnless(r_has("jsonlite"), "needs R with jsonlite")
    def test_an_unreadable_record_is_said_not_read_as_no_reasons(self):
        rec, err = self.record('{"excluded_runs": {"runs": [{"run": "MICH_03_rerun", "rea')
        self.assertEqual([x["reason"] for x in rec["runs"]], [None, None])
        self.assertIn(os.path.realpath(self.dec), rec["reasons_unreadable"])
        self.assertIn("could not be read (", rec["reasons_unreadable"])
        self.assertEqual(rec["reasons_from"], rec["reasons_unreadable"])
        self.assertIn("[run_de] WARNING: ", err)
        self.assertIn("could not be read", err)
        self.assertIn("2 run(s) left out of the DE are not used", err)
        self.assertTrue(all("reason NOT READ: the decisions record could not be read [confirm]" in m
                            for m in rec["methods"][1:]), rec["methods"])

    @unittest.skipUnless(r_has("jsonlite"), "needs R with jsonlite")
    def test_a_readable_record_still_gives_its_reasons(self):
        p = validate(self.meta, write_report(self.tmp, ANALYSED + RERUNS), *exclude_args(REASONS))
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        rec, err = self.record(None)
        self.assertEqual({x["run"]: x["reason"] for x in rec["runs"]}, REASONS)
        self.assertNotIn("reasons_unreadable", rec)
        self.assertNotIn("WARNING", err)

    @unittest.skipUnless(r_has("jsonlite"), "needs R with jsonlite")
    def test_no_record_and_an_odd_one_are_not_reasons_and_never_crash(self):
        rec, err = self.record(None)
        self.assertIn("--exclude was not used", rec["reasons_from"])
        self.assertNotIn("reasons_unreadable", rec)
        self.assertTrue(all("reason NOT RECORDED [confirm]" in m for m in rec["methods"][1:]))
        rec, err = self.record('{"excluded_runs": {"runs": ["MICH_03_rerun", 7]}}')
        self.assertEqual([x["reason"] for x in rec["runs"]], [None, None])
        self.assertNotIn("reasons_unreadable", rec)


if __name__ == "__main__":
    unittest.main()
