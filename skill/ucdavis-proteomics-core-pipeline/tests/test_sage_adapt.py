#!/usr/bin/env python3
"""
adapt_sage() turned Sage's lfq.parquet into the DE input with EVERY row: the decoys too, and
targets at any q_value. Sage's parquet writer (sage-cloudpath parquet.rs serialize_lfq, v0.14.7)
keeps the decoy MS1 peaks -- the same peptide in a +11.06 Da window, with the SAME `proteins`
string -- while its lfq.tsv writer (sage-cli output.rs write_lfq) drops them; and Sage counts an
MS1 peak as discovered only at q_value <= 0.05 (fdr.rs picked_precursor). The comment said
"Sage already FDR-filtered at write time". It had not.

Now a row reaches report.parquet only when is_decoy is false and q_value <= LFQ_Q_MAX, and what
was kept and dropped per file is recorded (sage_adapt.json, search_provenance.json).
"""
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

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False

# (peptide, proteins, is_decoy, q_value, {file: intensity}) -- Sage writes one row per
# precursor per file, and the decoy of a peptide carries the target's proteins string. Each
# target and its decoy share ONE q_value, as Sage stores it (fdr.rs picked_precursor keys the q
# map by precursor: 1,489 of 1,489 pairs identical in gabrig's HeL50 UnvPe lfq.parquet).
ROWS = [
    ("PEPTIDEAK", "sp|P1|A_HUMAN", False, 0.001, {"r1.mzML": 100.0, "r2.mzML": 110.0}),
    ("PEPTIDEAK", "sp|P1|A_HUMAN", True, 0.001, {"r1.mzML": 900.0, "r2.mzML": 950.0}),
    ("PEPTIDECK", "sp|P2|B_HUMAN", False, 0.05, {"r1.mzML": 50.0, "r2.mzML": 0.0}),    # at the line
    ("PEPTIDECK", "sp|P2|B_HUMAN", True, 0.05, {"r1.mzML": 40.0, "r2.mzML": 45.0}),
    ("PEPTIDEDK", "sp|P3|C_HUMAN", False, 0.128, {"r1.mzML": 70.0, "r2.mzML": 75.0}),  # fails
    ("PEPTIDEDK", "sp|P3|C_HUMAN", True, 0.128, {"r1.mzML": 71.0, "r2.mzML": 74.0}),
]


def write_lfq(path, rows=ROWS, drop=()):
    cols = {k: [] for k in ("peptide", "stripped_peptide", "charge", "proteins", "is_decoy",
                            "q_value", "filename", "intensity")}
    for pep, prot, dec, q, per in rows:
        for f, v in per.items():
            cols["peptide"].append(pep)
            cols["stripped_peptide"].append(pep)
            cols["charge"].append(None)
            cols["proteins"].append(prot)
            cols["is_decoy"].append(dec)
            cols["q_value"].append(q)
            cols["filename"].append(f)
            cols["intensity"].append(v)
    t = {k: v for k, v in cols.items() if k not in drop}
    if "charge" in t:
        t["charge"] = pa.array(t["charge"], pa.int32())
    for k in ("q_value", "intensity"):
        if k in t:
            t[k] = pa.array(t[k], pa.float32())
    pq.write_table(pa.table(t), path)


@unittest.skipUnless(HAVE_ARROW, "pyarrow is needed for the Sage fixtures")
class AdaptSageKeepsOnlyValidMs1Rows(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        import run_search
        self.rs = run_search

    def adapt(self):
        report = self.rs.adapt_sage(self.d)
        t = pq.read_table(report).to_pydict()
        return sorted(zip(t["Run"], t["Protein.Group"], t["PG.MaxLFQ"]))

    def test_exactly_the_target_rows_at_or_under_q_005(self):
        write_lfq(os.path.join(self.d, "lfq.parquet"))
        got = self.adapt()
        # Protein.Group as DIA-NN writes it -- bare accessions (protein_ids.group_accessions) --
        # so the contaminant filter's "starts with Cont_" rule can see Sage's contaminants
        self.assertEqual(got, sorted([
            ("r1", "P1", 100.0), ("r2", "P1", 110.0),
            ("r1", "P2", 50.0), ("r2", "P2", 0.0)]))
        names = pq.read_table(os.path.join(self.d, "report.parquet")).column("Protein.Names")
        self.assertEqual(set(names.to_pylist()), {"A_HUMAN", "B_HUMAN"})
        # none of the decoy intensities, none of the failing target's
        vals = {v for _, _, v in got}
        self.assertFalse(vals & {900.0, 950.0, 40.0, 45.0, 70.0, 75.0, 71.0, 74.0})

    def test_per_file_counts_are_recorded_and_merged_into_provenance(self):
        write_lfq(os.path.join(self.d, "lfq.parquet"))
        with open(os.path.join(self.d, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "sage"}, fh)
        self.adapt()
        with open(os.path.join(self.d, "sage_adapt.json")) as fh:
            rec = json.load(fh)
        want = {"rows": 6, "decoy_dropped": 3, "q_value_dropped": 1, "kept": 2}
        self.assertEqual(rec["per_file"], {"r1": want, "r2": want})
        self.assertEqual(rec["totals"], {k: 2 * v for k, v in want.items()})
        self.assertEqual(rec["q_value_max"], 0.05)
        self.assertIn("write_lfq", rec["decoy_rule_source"])
        self.assertIn("picked_precursor", rec["q_value_source"])
        with open(os.path.join(self.d, "search_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov["engine"], "sage")
        self.assertEqual(prov["sage_adapt"]["totals"]["decoy_dropped"], 6)

    def test_sage_s_own_count_is_recorded_beside_what_was_kept(self):
        """The stored q is shared with the decoy, so the kept precursors are a conservative
        subset of Sage's logged count -- both are on the record, never one passed off as the
        other."""
        write_lfq(os.path.join(self.d, "lfq.parquet"))
        with open(os.path.join(self.d, "sage.log"), "w") as fh:
            fh.write("[INFO sage] discovered 3 target MS1 peaks at 5% FDR\n")
        self.adapt()
        with open(os.path.join(self.d, "sage_adapt.json")) as fh:
            rec = json.load(fh)
        self.assertEqual(rec["target_precursors_kept"], 2)          # PEPTIDEAK, PEPTIDECK
        self.assertEqual(rec["sage_logged_target_ms1_peaks_5pct"], 3)
        self.assertIn("conservative subset", rec["kept_rule"])
        self.assertIn("shared by the target and its decoy", rec["q_value_source"])

    def test_the_threshold_is_the_one_definition(self):
        import sage_lfq_check
        write_lfq(os.path.join(self.d, "lfq.parquet"))
        self.addCleanup(setattr, sage_lfq_check, "LFQ_Q_MAX", sage_lfq_check.LFQ_Q_MAX)
        sage_lfq_check.LFQ_Q_MAX = 0.2
        self.assertEqual(len(self.adapt()), 6)       # the q 0.128 target now passes; decoys never

    def test_a_file_without_decoy_labels_is_refused(self):
        write_lfq(os.path.join(self.d, "lfq.parquet"), drop=("is_decoy",))
        with self.assertRaises(SystemExit) as cm:
            self.rs.adapt_sage(self.d)
        self.assertIn("cannot be told from the targets", str(cm.exception))
        self.assertFalse(os.path.exists(os.path.join(self.d, "report.parquet")))

    def test_nothing_passing_is_said_not_passed_off(self):
        """gabrig's +7 ppm cohort: every target at q 0.128. The old adapter still handed the DE
        a full matrix of decoy and failing-target intensities."""
        write_lfq(os.path.join(self.d, "lfq.parquet"),
                  rows=[r for r in ROWS if r[3] == 0.128])
        r = subprocess.run([sys.executable, "-c",
                            "import sys; sys.path.insert(0, sys.argv[1]); import run_search; "
                            "run_search.adapt_sage(sys.argv[2])", SCRIPTS, self.d],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("no Sage MS1 peak passed q_value <= 0.05", r.stderr)
        self.assertEqual(pq.read_table(os.path.join(self.d, "report.parquet")).num_rows, 0)


if __name__ == "__main__":
    unittest.main()
