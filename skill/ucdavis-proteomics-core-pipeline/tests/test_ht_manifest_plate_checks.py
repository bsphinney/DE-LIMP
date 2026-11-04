"""
Three plate checks in `ht_manifest.py fetch` (I41, 2.14.2; Brett's Q7: the checks only, no queue
building), each a hard stop on FAIL like the other gates:

* `run_sizes` -- a blank of the plate (STAN's class, or STAN's blank words in its name), larger
  than BLANK_FRACTION of the plate's median sample run, is a FAIL: the run names no longer match
  what is in the runs (on two earlier HT submissions, 94 of 96 runs were named for the wrong well and passed
  every gate; runs named Blank_S3-A6 held 1.36 GB). STAN leaves blanks out of its list with the
  default --include samples, so they are looked for beside the plate's runs. A sample run under
  SMALL_RUN_FRACTION of the median is a WARN (a failed injection?), and so is a plate with no
  blank to check.
* `sample_list_copied` -- the submission's CoreOmics sample list is the same ordered list as
  another submission's (8 or more names): a list copied from another plate.
* `sample_table` -- the runs' HyStar sample table (SampleInfo.xml SampleTablePath): more than one
  table, or a table naming ANOTHER PROT number. Real table names are typed by staff and carry no
  PROT number; such a name is OK, and named in the gate (controller ruling 34, 2026-10-31).

SampleInfo.xml is built here from the layout of a real one (UTF-16 little-endian with a BOM, no line
breaks; <SampleTable> holding <SampleTableHeader ... SampleTablePath=...>, <Sample ...>,
<AnalysisHeader .../> and <SampleTableProperties>), with made-up names. Run sizes come from
truncated files, all small: the limits are relative to the plate, and the Mac's DataArchive disk is
HFS+, which has no sparse files (a 1.2 GB truncate allocates 1.2 GB there).

Stdlib only; no STAN, no network: the fake `stan` of test_ht_manifest_gates.
"""
import argparse
import codecs
import json
import os
import sys
import tempfile
import unittest
from unittest import mock
from xml.sax.saxutils import quoteattr

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import ht_manifest  # noqa: E402
from test_ht_manifest_gates import _plate_run, _run, _well  # noqa: E402

SAMPLE = 200_000            # a sample run's analysis.tdf_bin, scaled down: only ratios matter
TDF = 10_000
BLANK = 1_000               # a real blank: far smaller than a sample run
TABLE = "Plate #3 Example words October 27th 2026"     # as staff type it: no PROT number
PREFIX = "20261023_PROT_9001_100spd"


def sample_info_xml(table, data_path, position):
    """A per-run SampleInfo.xml in the real layout, encoded as the instrument writes it."""
    def attrs(d):
        return " ".join(f"{k}={quoteattr(str(v))}" for k, v in d.items())
    header = {"Version": "2.0", "CreationDate": "2026-10-23T09:12:00", "Creator": "staffexample",
              "Description": "", "NumOfSamples": "1", "HyStarVersion": "6.2.1.13",
              "SampleTablePath": table, "BatchID": "1"}
    sample = {"analysisID": "1", "Line": "1", "Position": position, "SampleID": "XY1",
              "DataPath": data_path, "SampleType": "Sample", "ResultDatafile": data_path,
              "Method": "D:\\Methods\\100spd.m", "MS_Method": "D:\\Methods\\100spd.m",
              "Description": "", "Submitter_Name": "Example Submitter",
              "Submitter_Lab": "Example Lab"}
    analysis = {"CreationDateTime": "2026-10-23T09:12:00", "FileName": data_path,
                "SampleID": "XY1", "Method": "100spd.m"}
    props = "".join(f"<Property {attrs({'PropertyNumber': i, 'Name': f'Property{i}', 'Value': ''})}/>"
                    for i in range(1, 9))
    xml = ('<?xml version="1.0" encoding="utf-16"?>'
           f"<SampleTable><SampleTableHeader {attrs(header)}/><Sample {attrs(sample)}/>"
           f"<AnalysisHeader {attrs(analysis)}/><SampleTableProperties>{props}"
           "</SampleTableProperties></SampleTable>")
    return codecs.BOM_UTF16_LE + xml.encode("utf-16-le")


def d_run(tmp, i, size=SAMPLE, sample=None, table=None, prefix=PREFIX, inj=25000, plate="S3",
          name=None):
    """A timsTOF run folder of the plate: analysis.tdf_bin of `size` bytes and analysis.tdf, and
    with `table` its SampleInfo.xml."""
    name = name or f"{prefix}_{sample or f'XY{i + 1}'}_{plate}-{_well(i)}_1_{inj + i}.d"
    d = os.path.join(tmp, name)
    os.makedirs(d)
    for fn, n in (("analysis.tdf_bin", size), ("analysis.tdf", TDF)):
        p = os.path.join(d, fn)
        open(p, "w").close()
        os.truncate(p, n)
    if table is not None:
        with open(os.path.join(d, "SampleInfo.xml"), "wb") as fh:
            fh.write(sample_info_xml(table, "D:\\Data\\Oct26\\" + name, f"{plate}-{_well(i)}"))
    return d


def manifest(files, entries=None):
    m = {"submission": "9001", "include": "samples", "files": files, "n_files": len(files),
         "plates": ["S3"], "counts": {"sample": len(files)}, "n_needs_rerun": 0,
         "missing_paths": []}
    if entries is not None:
        m["entries"] = entries
    return m


def record(folder, names, internal_id):
    """core_submission.py fetch's submission_summary.json, as much of it as the checks read."""
    os.makedirs(folder, exist_ok=True)
    p = os.path.join(folder, "submission_summary.json")
    with open(p, "w") as fh:
        json.dump({"internal_id": internal_id, "id": "0123456789ab",
                   "samples": [{"unique_id": f"U{i + 1:03d}", "sample_name": n,
                                "condition_name": "b" if i % 2 else "a"}
                               for i, n in enumerate(names)]}, fh)
    return p


def gates(tmp):
    with open(os.path.join(tmp, "ht_manifest.json")) as fh:
        return {g["gate"]: g for g in json.load(fh)["gates"]}


class RunSizes(unittest.TestCase):
    def test_blank_the_size_of_a_sample_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i) for i in range(12)]
            big = d_run(tmp, 12, sample="Blank1")
            small = d_run(tmp, 13, size=BLANK, sample="Blank2")
            r = _run(tmp, manifest(runs + [big, small]))
            self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
            g = gates(tmp)["run_sizes"]
            self.assertEqual(g["status"], "FAIL")
            self.assertIn(os.path.basename(big), g["detail"] + " ".join(g["listed"]))
            self.assertNotIn(os.path.basename(small), " ".join(g["listed"]))
            self.assertIn(os.path.basename(big), r.stdout)
            self.assertIn("HARD GATE FAILED", r.stderr)

    def test_a_blank_stan_left_out_is_found_beside_the_runs(self):
        """STAN classes a name with its blank words as a blank and, with --include samples, leaves
        it out of `files` and `entries`: one plate's sample-sized runs named Blank_S3-A6 were never
        measured. They are looked for in the plate runs' folders (one listing per folder)."""
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i) for i in range(12)]
            big = d_run(tmp, 12, sample="Blank")
            small = d_run(tmp, 13, size=BLANK, sample="Wash")
            other = d_run(tmp, 14, prefix="20261023_PROT_9002_100spd", sample="Blank")
            r = _run(tmp, manifest(runs, [{"raw_path": f, "class": "sample"} for f in runs]))
            self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
            g = gates(tmp)["run_sizes"]
            self.assertEqual(g["status"], "FAIL")
            listed = " ".join(g["listed"])
            self.assertIn(os.path.basename(big), listed)
            self.assertNotIn(os.path.basename(small), listed)
            self.assertNotIn(os.path.basename(other), listed)       # another submission's blank
            self.assertEqual(g["n_blanks"], 2)
            self.assertNotIn(big, open(os.path.join(tmp, "files.txt")).read())

    def test_a_blank_by_stans_class_counts_too(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i) for i in range(4)]
            entries = [{"raw_path": f, "class": "sample"} for f in runs]
            entries[2]["class"] = "blank"
            r = _run(tmp, manifest(runs, entries))
            self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
            g = gates(tmp)["run_sizes"]
            self.assertEqual(g["status"], "FAIL")
            self.assertIn(os.path.basename(runs[2]), " ".join(g["listed"]))

    def test_stans_blank_words(self):
        for name in ("x_Blank_S3-A1", "x_wash2_S3-A1", "x_BLK_S3-A1", "x_blnk_S3-A1"):
            self.assertTrue(ht_manifest.BLANK_NAME.search(name), name)
        self.assertIsNone(ht_manifest.BLANK_NAME.search("x_Liver3_S3-A1"))

    def test_the_limit_is_relative_to_the_plate(self):
        """A fixed size would fail real blanks at a lower throughput, where every run is bigger."""
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i) for i in range(6)]
            half = d_run(tmp, 6, size=SAMPLE // 2 + TDF, sample="Blank")   # just over half
            plate = ht_manifest.submission_plate("9001")
            g = ht_manifest.run_size_gate(runs, {}, plate)
            self.assertEqual(g["status"], "FAIL", g)
            self.assertIn(os.path.basename(half), " ".join(g["listed"]))
            with mock.patch.object(ht_manifest, "BLANK_FRACTION", 0.9):
                self.assertEqual(ht_manifest.run_size_gate(runs, {}, plate)["status"], "PASS")

    def test_small_run_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i) for i in range(12)] + [d_run(tmp, 13, size=BLANK, sample="Blank")]
            tiny = d_run(tmp, 12, size=SAMPLE // 10)
            r = _run(tmp, manifest(runs + [tiny]))
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            g = gates(tmp)["run_sizes"]
            self.assertEqual(g["status"], "WARN")
            self.assertEqual(len(g["small"]), 1)
            self.assertIn(os.path.basename(tiny), " ".join(g["listed"]))
            self.assertIn("failed injection", g["detail"])
            self.assertIn(os.path.basename(tiny), r.stdout)

    def test_no_blank_to_check_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i) for i in range(12)]
            self.assertEqual(_run(tmp, manifest(runs)).returncode, 0)
            g = gates(tmp)["run_sizes"]
            self.assertEqual(g["status"], "WARN")
            self.assertIn("no blank", g["detail"])

    def test_a_plate_of_even_runs_and_small_blanks_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i) for i in range(12)]
            d_run(tmp, 12, size=BLANK, sample="Blank1")                 # beside them, not listed
            self.assertEqual(_run(tmp, manifest(runs)).returncode, 0)
            g = gates(tmp)["run_sizes"]
            self.assertEqual((g["status"], g["n_blanks"]), ("PASS", 1), g)


class CopiedSampleList(unittest.TestCase):
    NAMES = [f"XY{i + 1}" for i in range(96)]        # the plate runs are named XY<n>

    def plate(self, tmp):
        return manifest([_plate_run(tmp, i, prefix=PREFIX, plate="S3", inj=25000)
                         for i in range(96)])

    def test_sample_list_identical_to_another_submission_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            core = os.path.join(tmp, "core")
            ours = record(os.path.join(core, "PROT_9001"), self.NAMES, "PROT_9001")
            record(os.path.join(core, "PROT_9793"), self.NAMES, "PROT_9793")
            r = _run(tmp, self.plate(tmp), ["--summary", ours])
            self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
            g = gates(tmp)["sample_list_copied"]
            self.assertEqual(g["status"], "FAIL")
            self.assertIn("PROT_9793", g["detail"])
            self.assertEqual(g["copied_from"], ["PROT_9793"])
            self.assertIn("PROT_9793", r.stdout)

    def test_another_order_or_a_short_list_is_no_copy(self):
        with tempfile.TemporaryDirectory() as tmp:
            core = os.path.join(tmp, "core")
            ours = record(os.path.join(core, "PROT_9001"), self.NAMES, "PROT_9001")
            record(os.path.join(core, "PROT_9793"), list(reversed(self.NAMES)), "PROT_9793")
            record(os.path.join(core, "PROT_9794"), self.NAMES[:40], "PROT_9794")
            record(os.path.join(core, "PROT_9001_copy"), self.NAMES, "PROT_9001")   # our own
            r = _run(tmp, self.plate(tmp), ["--summary", ours])
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            g = gates(tmp)["sample_list_copied"]
            self.assertEqual(g["status"], "PASS", g)
            self.assertEqual(g["n_compared"], 2)

    def test_fewer_than_eight_samples_are_not_compared(self):
        with tempfile.TemporaryDirectory() as tmp:
            core = os.path.join(tmp, "core")
            names = self.NAMES[:6]
            ours = record(os.path.join(core, "PROT_9001"), names, "PROT_9001")
            record(os.path.join(core, "PROT_9793"), names, "PROT_9793")
            m = manifest([_plate_run(tmp, i, prefix=PREFIX, plate="S3", inj=25000) for i in range(6)])
            r = _run(tmp, m, ["--summary", ours])
            self.assertEqual(gates(tmp)["sample_list_copied"]["status"], "PASS", r.stdout)

    def test_summaries_dir_names_where_the_other_records_are(self):
        with tempfile.TemporaryDirectory() as tmp:
            ours = record(os.path.join(tmp, "here"), self.NAMES, "PROT_9001")
            elsewhere = os.path.join(tmp, "elsewhere")
            record(os.path.join(elsewhere, "PROT_9793"), self.NAMES, "PROT_9793")
            r = _run(tmp, self.plate(tmp), ["--summary", ours, "--summaries-dir", elsewhere])
            self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
            self.assertIn("PROT_9793", gates(tmp)["sample_list_copied"]["detail"])

    def test_without_the_record_it_is_not_checked_and_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(_run(tmp, self.plate(tmp)).returncode, 0)
            g = gates(tmp)["sample_list_copied"]
            self.assertEqual(g["status"], "WARN")
            self.assertIn("not checked", g["detail"])

    # ---- review of 9d2559a5: false copies, and one bad record must not stop the fetch --------
    def gate(self, ours, summaries_dir=None):
        return ht_manifest.copied_list_gate(argparse.Namespace(summary=ours,
                                                               summaries_dir=summaries_dir))

    def test_generic_names_are_not_compared(self):
        """Two unrelated submissions both named their samples 1..96 (or A1..H12, Sample 1..):
        the same list by convention, not a copy."""
        for names in ([str(i + 1) for i in range(96)],
                      [f"{r}{c}" for r in "ABCDEFGH" for c in range(1, 13)],
                      [f"Sample {i + 1}" for i in range(12)]):
            with self.subTest(names[0]), tempfile.TemporaryDirectory() as tmp:
                core = os.path.join(tmp, "core")
                ours = record(os.path.join(core, "PROT_9001"), names, "PROT_9001")
                record(os.path.join(core, "PROT_9793"), names, "PROT_9793")
                g = self.gate(ours)
                self.assertEqual(g["status"], "INFO", g)
                self.assertIn("not compared: the names are not distinctive", g["detail"])

    def test_a_list_with_empty_names_is_not_compared(self):
        with tempfile.TemporaryDirectory() as tmp:
            core = os.path.join(tmp, "core")
            ours = record(os.path.join(core, "PROT_9001"), [""] * 12, "PROT_9001")
            record(os.path.join(core, "PROT_9793"), [""] * 12, "PROT_9793")
            g = self.gate(ours)
            self.assertEqual(g["status"], "INFO", g)
            self.assertIn("12 sample(s) have no name", g["detail"])

    def test_a_partly_generic_list_is_still_compared(self):
        """(pin) Distinctive names beside numbers: still a copy when the whole list matches."""
        names = [str(i + 1) for i in range(10)] + ["Liver_ctrl", "Liver_treat"]
        with tempfile.TemporaryDirectory() as tmp:
            core = os.path.join(tmp, "core")
            ours = record(os.path.join(core, "PROT_9001"), names, "PROT_9001")
            record(os.path.join(core, "PROT_9793"), names, "PROT_9793")
            self.assertEqual(self.gate(ours)["status"], "FAIL")

    def test_a_malformed_record_is_named_not_a_crash(self):
        with tempfile.TemporaryDirectory() as tmp:
            core = os.path.join(tmp, "core")
            ours = record(os.path.join(core, "PROT_9001"), self.NAMES, "PROT_9001")
            bad = os.path.join(core, "PROT_9793", "submission_summary.json")
            os.makedirs(os.path.dirname(bad))
            with open(bad, "w") as fh:
                json.dump({"internal_id": "PROT_9793", "samples": 5}, fh)
            r = _run(tmp, self.plate(tmp), ["--summary", ours])
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            g = gates(tmp)["sample_list_copied"]
            self.assertEqual(g["status"], "INFO")
            self.assertIn(bad, " ".join(g["listed"]))

    def test_our_own_malformed_record_is_not_a_crash(self):
        with tempfile.TemporaryDirectory() as tmp:
            ours = os.path.join(tmp, "core", "PROT_9001", "submission_summary.json")
            os.makedirs(os.path.dirname(ours))
            with open(ours, "w") as fh:
                json.dump({"internal_id": "PROT_9001", "samples": 5}, fh)
            self.assertEqual(self.gate(ours)["status"], "WARN")


OTHER = "Plate #4 Other words October 28th 2026"
RERUN = "20261030_PROT_9001_rerun_XY4_S1-A1_1_25200.d"     # a re-injection, on another tray


class SampleTable(unittest.TestCase):
    """One table per PLATE (controller ruling 34b, 2026-10-31: HIVE's Sep/Oct 2026 HT plates,
    four submissions, have 96 runs and exactly one SampleTablePath each, so a two-tray
    submission is most likely two tables). A re-injection, or a run whose plate is unknown, is
    checked only for a table naming another submission."""

    def test_runs_from_two_sample_tables_fail(self):
        """One plate (S3), two tables."""
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i, table=TABLE if i < 6 else OTHER) for i in range(12)]
            r = _run(tmp, manifest(runs))
            self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
            g = gates(tmp)["sample_table"]
            self.assertEqual(g["status"], "FAIL")
            self.assertEqual(sorted(g["tables"]), sorted([TABLE, OTHER]))
            self.assertIn(TABLE, g["detail"])
            self.assertIn(OTHER, g["detail"])
            self.assertIn("plate S3", g["detail"])
            self.assertIn("core_submission.py locate --files-from", g["detail"])

    def test_two_plates_with_a_table_each_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = ([d_run(tmp, i, table=TABLE) for i in range(6)]
                    + [d_run(tmp, i, table=OTHER, plate="S4", inj=25100) for i in range(6)])
            r = _run(tmp, manifest(runs))
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            g = gates(tmp)["sample_table"]
            self.assertEqual(g["status"], "PASS", g)
            self.assertEqual(g["plate_tables"], {"S3": [TABLE], "S4": [OTHER]})
            self.assertIn(f"S3: {TABLE}", r.stdout)
            self.assertIn(f"S4: {OTHER}", r.stdout)

    def test_the_plate_is_stans_when_stan_gives_one(self):
        """STAN's entry `plate` before the name's _S<n>- tail."""
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i, table=TABLE if i < 3 else OTHER) for i in range(6)]
            entries = [{"raw_path": f, "class": "sample", "plate": "S3" if i < 3 else "S4"}
                       for i, f in enumerate(runs)]
            self.assertEqual(_run(tmp, manifest(runs, entries)).returncode, 0)
            self.assertEqual(gates(tmp)["sample_table"]["status"], "PASS")

    def test_a_reinjection_queued_from_its_own_table_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i, table=TABLE) for i in range(6)]
            runs.append(d_run(tmp, 0, table="Reruns October 30th 2026", name=RERUN))
            r = _run(tmp, manifest(runs))
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertEqual(gates(tmp)["sample_table"]["status"], "PASS")

    def test_a_reinjection_from_another_submissions_table_fails(self):
        """(pin) A re-injection is still checked for a table naming another submission."""
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i, table=TABLE) for i in range(6)]
            runs.append(d_run(tmp, 0, table="PROT_9838 reruns", name=RERUN))
            r = _run(tmp, manifest(runs))
            self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
            g = gates(tmp)["sample_table"]
            self.assertEqual(g["status"], "FAIL")
            self.assertIn("PROT_9838", g["detail"])

    def test_a_table_naming_this_and_another_submission_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i, table="PROT_9001 and PROT_9838 plate 3") for i in range(6)]
            r = _run(tmp, manifest(runs))
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            g = gates(tmp)["sample_table"]
            self.assertEqual(g["status"], "WARN")
            self.assertIn("PROT_9838", g["detail"])

    def test_table_named_for_another_submission_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            foreign = "PROT_9838 Plate #3 Example words"
            runs = [d_run(tmp, i, table=foreign) for i in range(6)]
            r = _run(tmp, manifest(runs))
            self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
            g = gates(tmp)["sample_table"]
            self.assertEqual(g["status"], "FAIL")
            self.assertIn("PROT_9838", g["detail"])

    def test_a_table_with_no_submission_number_passes_and_is_named(self):
        """Ruling 34: real HyStar table names carry no PROT number."""
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i, table=TABLE) for i in range(6)]
            r = _run(tmp, manifest(runs))
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            g = gates(tmp)["sample_table"]
            self.assertEqual(g["status"], "PASS")
            self.assertIn(TABLE, g["detail"])
            self.assertIn(f"sample table: {TABLE}", r.stdout)

    def test_a_table_naming_this_submission_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i, table="PROT_9001 plate 1") for i in range(6)]
            self.assertEqual(_run(tmp, manifest(runs)).returncode, 0)
            self.assertEqual(gates(tmp)["sample_table"]["status"], "PASS")

    def test_without_sample_info_it_is_not_checked_and_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i) for i in range(6)]
            self.assertEqual(_run(tmp, manifest(runs)).returncode, 0)
            g = gates(tmp)["sample_table"]
            self.assertEqual(g["status"], "WARN")
            self.assertIn("not checked", g["detail"])

    def test_an_unreadable_sample_info_is_named(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = [d_run(tmp, i, table=TABLE) for i in range(6)]
            with open(os.path.join(runs[3], "SampleInfo.xml"), "wb") as fh:
                fh.write(b"\xff\xfe<Sam")
            self.assertEqual(_run(tmp, manifest(runs)).returncode, 0)
            g = gates(tmp)["sample_table"]
            self.assertEqual(g["status"], "WARN")
            self.assertIn(os.path.basename(runs[3]), " ".join(g["listed"]))

    def test_the_table_name_is_read_as_the_instrument_writes_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = d_run(tmp, 0, table=TABLE)
            with open(os.path.join(run, "SampleInfo.xml"), "rb") as fh:
                self.assertTrue(fh.read(2) == codecs.BOM_UTF16_LE)
            self.assertEqual(ht_manifest.sample_table_path(run), TABLE)


class TheExitContract(unittest.TestCase):
    """The docstring's contract: 2 = HARD GATE FAILED, do not search -- for each new check."""

    def test_fetch_exit_code_is_2_on_any_of_them(self):
        names = [f"XY{i + 1}" for i in range(12)]

        def blank(tmp):
            d_run(tmp, 4, sample="Blank1")                    # beside the runs, as STAN leaves it
            return manifest([d_run(tmp, i) for i in range(4)]), []

        def copied(tmp):
            core = os.path.join(tmp, "core")
            ours = record(os.path.join(core, "PROT_9001"), names, "PROT_9001")
            record(os.path.join(core, "PROT_9793"), names, "PROT_9793")
            return manifest([_plate_run(tmp, i, prefix=PREFIX, plate="S3", inj=25000)
                             for i in range(12)]), ["--summary", ours]

        def two_tables(tmp):
            return manifest([d_run(tmp, i, table=TABLE if i else "Plate #9 Elsewhere")
                             for i in range(4)]), []

        def foreign_table(tmp):
            return manifest([d_run(tmp, i, table="PROT_9838 plate") for i in range(4)]), []

        for why, build in (("run_sizes", blank), ("sample_list_copied", copied),
                           ("sample_table", two_tables), ("sample_table", foreign_table)):
            with self.subTest(build.__name__), tempfile.TemporaryDirectory() as tmp:
                m, extra = build(tmp)
                r = _run(tmp, m, extra)
                self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
                self.assertEqual(gates(tmp)[why]["status"], "FAIL")
        for g in ("run_sizes", "sample_list_copied", "sample_table", "--summaries-dir"):
            self.assertIn(g, ht_manifest.__doc__)
        self.assertNotIn("400 MB", ht_manifest.__doc__)

    def test_the_reference_describes_them(self):
        ref = os.path.join(os.path.dirname(SCRIPTS), "references", "ht-submissions.md")
        with open(ref, encoding="utf-8") as fh:
            doc = fh.read()
        for g in ("`run_sizes`", "`sample_list_copied`", "`sample_table`", "--summaries-dir"):
            self.assertTrue(g in doc, f"references/ht-submissions.md does not describe {g}")


if __name__ == "__main__":
    unittest.main()
