#!/usr/bin/env python3
"""The search's quantification record has an engine-neutral name, and every reader takes the old
one too (2.15 F4; LEDGER: "the record file name sage_lfq_check.json is reused for FragPipe; an
engine-neutral rename is a 2.15 item").

  * sage_lfq_check.RECORD_NAMES: quantification_check.json first, then 2.14.x's
    sage_lfq_check.json; load() / record_path() read either, the newer first;
  * run_de.R gates on either name (its list is tied to RECORD_NAMES here) and refuses a DE beside
    an identifications-only record of either name;
  * a delivered copy of either names the role, never the login (staff.ACCEPTED_BY_RECORDS);
  * OUTPUT_FILES.md describes both;
  * an analysis delivery carries a FragPipe search's combined_*.tsv and fragpipe_quantity.json,
    SKIPPED when absent (T5).
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import sage_lfq_check as slc                    # noqa: E402
import staff                                    # noqa: E402
from test_run_de_contaminants import r_has      # noqa: E402

NEW, OLD = "quantification_check.json", "sage_lfq_check.json"
FRAGPIPE_TABLES = ("combined_protein.tsv", "combined_peptide.tsv",
                   "combined_modified_peptide.tsv", "combined_ion.tsv", "fragpipe_quantity.json")


def write(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(obj if isinstance(obj, str) else json.dumps(obj))
    return path


class Names(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def test_the_neutral_name_first_then_the_old(self):
        self.assertEqual(slc.RECORD_NAMES, (NEW, OLD))
        self.assertIn(slc.RECORD, slc.RECORD_NAMES)

    def test_either_name_is_read_the_newer_first(self):
        self.assertIsNone(slc.load(self.d))
        write(os.path.join(self.d, OLD), {"status": "not_applicable"})
        self.assertEqual(slc.record_path(self.d), os.path.join(self.d, OLD))
        self.assertFalse(slc.search_quantified(self.d))
        write(os.path.join(self.d, NEW), {"status": "ok"})
        self.assertEqual(slc.record_path(self.d), os.path.join(self.d, NEW))
        self.assertEqual(slc.load(self.d)["status"], "ok")

    def test_the_gate_answers_for_either_name(self):
        for name in (NEW, OLD):
            d = tempfile.mkdtemp(dir=self.d)
            write(os.path.join(d, name), {"status": "not_applicable",
                                          "de_refusal": "REFUSED: no quantities here"})
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sage_lfq_check.py"), "gate",
                                "--out", d], capture_output=True, text=True, timeout=60)
            j = json.loads(r.stdout)
            self.assertEqual((j["checked"], j["refused"]), (True, True), name)

    def test_run_de_reads_the_same_names(self):
        with open(os.path.join(SCRIPTS, "run_de.R")) as fh:
            text = fh.read()
        m = re.search(r'QUANTIFICATION_RECORDS <- c\(([^)]*)\)', text)
        self.assertIsNotNone(m)
        self.assertEqual(tuple(re.findall(r'"([^"]+)"', m.group(1))), slc.RECORD_NAMES)
        self.assertNotIn('file.path(search_dir, "sage_lfq_check.json")', text)

    def test_the_watcher_reads_the_same_names_newer_first(self):
        with open(os.path.join(SCRIPTS, "watch_run.sh")) as fh:
            line = next(ln for ln in fh if ln.startswith('[ -n "$OUTDIR" ] && lfq_rec='))
        self.assertEqual(tuple(re.findall(r'\$OUTDIR/([a-z_]+\.json)', line)), slc.RECORD_NAMES)

    def test_a_delivered_copy_of_either_names_the_role(self):
        for name in slc.RECORD_NAMES:
            self.assertEqual(staff.ACCEPTED_BY_RECORDS[name], ("accepted",))
            p = write(os.path.join(self.d, name),
                      {"status": "warn", "accepted": {"by": "someone123", "reason": "r"}})
            copy = staff.client_copy(p)
            self.assertIsNotNone(copy, name)
            self.assertNotIn("someone123", copy)

    def test_the_catalog_describes_both(self):
        import make_report
        for name in slc.RECORD_NAMES + FRAGPIPE_TABLES + (
                "directlfq_record.json", "combined_ion.tsv.protein_intensities.tsv",
                "combined_ion.tsv.no_norm.protein_intensities.tsv"):
            cat, desc = make_report.describe(name)
            self.assertNotIn("unrecognized", desc, name)
            self.assertEqual(cat, "Search output", name)
        self.assertNotIn("MaxLFQ table", make_report.describe("combined_protein.tsv")[1])


@unittest.skipUnless(r_has("jsonlite"), "needs R with jsonlite")
class RunDeGatesOnEitherName(unittest.TestCase):
    def test_a_de_beside_an_ids_only_record_is_refused_under_either_name(self):
        for name in (NEW, OLD):
            d = tempfile.mkdtemp()
            self.addCleanup(shutil.rmtree, d, True)
            out = os.path.join(d, "search")
            write(os.path.join(out, name), {
                "status": "not_applicable",
                "de_refusal": "REFUSED: this FragPipe search ran no quantification (test)"})
            write(os.path.join(out, "search_provenance.json"), {"engine": "fragpipe"})
            report = os.path.join(out, "report.parquet")
            open(report, "w").close()                  # the gate stops before reading it
            cond = write(os.path.join(d, "conditions.csv"),
                         "File.Name,Group\nS1,A\nS2,A\nS3,B\nS4,B\n")
            env = dict(os.environ)
            env["PATH"] = os.pathsep.join([os.path.dirname(sys.executable),
                                           os.path.dirname(shutil.which("Rscript")),
                                           "/usr/bin", "/bin"])
            r = subprocess.run(["Rscript", os.path.join(SCRIPTS, "run_de.R"), "--input", report,
                                "--metadata", cond, "--method", "maxlfq",
                                "--outdir", os.path.join(d, "de")],
                               capture_output=True, text=True, timeout=300, cwd=d, env=env)
            self.assertNotEqual(r.returncode, 0, name)
            self.assertIn("ran no quantification (test)", r.stderr, name)


class FragPipeTablesAreDelivered(unittest.TestCase):
    def test_the_delivery_plan_takes_them(self):
        import core_submission
        with tempfile.TemporaryDirectory() as d:
            for n in FRAGPIPE_TABLES:
                write(os.path.join(d, "search", n), "x\n")
            items = {i["name"]: i for i in core_submission.plan_delivery(d)}
        for n in FRAGPIPE_TABLES:
            self.assertEqual(items["search/" + n]["status"], "OK", items.get("search/" + n))

    def test_a_search_without_them_skips_them(self):
        import core_submission
        with tempfile.TemporaryDirectory() as d:
            items = {i["name"]: i for i in core_submission.plan_delivery(d)}
        for n in FRAGPIPE_TABLES:
            self.assertEqual(items["search/" + n]["status"], "SKIPPED")
            self.assertFalse(items["search/" + n]["required"])

    def test_the_fragpipe_logs_and_manifests_stay_behind(self):
        import core_submission
        for n in ("fragpipe.workflow", "fragpipe.fp-manifest", "fragpipe-files.fp-manifest",
                  "experiment_annotation.tsv", "log_2026-11-01_20-59-15.txt"):
            self.assertNotIn(n, core_submission.SEARCH_FILES)


if __name__ == "__main__":
    unittest.main()
