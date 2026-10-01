#!/usr/bin/env python3
"""
FragPipe vs DIA-NN on the same raw files shared 0 proteins in compare_searches.py.

run_search.adapt_fragpipe_dda() took FragPipe's `Protein` column (sp|P12345|ALBU_HUMAN) as the
Protein.Group, where DIA-NN reports the accession (P12345), and compare_searches.py compared the
two strings as they were. Now the adapter takes `Protein ID` (keeping FragPipe's own
`contam_` tag, which `Protein ID` drops), and compare_searches.py compares every report on its
accession (protein_ids.normalize_protein_id), so a report adapted before the fix still meets
DIA-NN's.

The FragPipe fixture has the column layout of a real IonQuant combined_protein.tsv (FragPipe's
Protein, Protein ID, Entry Name, Gene, ..., "<sample> MaxLFQ Intensity"), including both
contaminant spellings seen in Core runs: FragPipe's `contam_sp|...` and the skill's
`sp|Cont_...`.
"""
import csv
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
import protein_ids as pid     # noqa: E402

try:
    import pyarrow.parquet as pq      # noqa: F401
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False

RUNS = ("S1", "S2", "S3")
# (Protein, Protein ID, Entry Name, Gene) -- as FragPipe writes them
FRAGPIPE = [
    ("sp|P02768|ALBU_HUMAN", "P02768", "ALBU_HUMAN", "ALB"),
    ("sp|P60709|ACTB_HUMAN", "P60709", "ACTB_HUMAN", "ACTB"),
    ("tr|A0A075B6S2|A0A075B6S2_HUMAN", "A0A075B6S2", "A0A075B6S2_HUMAN", "IGKV2D-29"),
    ("sp|P04637-2|P53_HUMAN", "P04637-2", "P53_HUMAN", "TP53"),
    ("sp|Q99999|ONLYFP_HUMAN", "Q99999", "ONLYFP_HUMAN", "ONLYFP"),
    ("contam_sp|P00167|CYB5_HUMAN", "P00167", "CYB5_HUMAN", "CYB5A"),
    ("sp|Cont_P00761|TRYP_PIG", "Cont_P00761", "TRYP_PIG", ""),
]
# DIA-NN's Protein.Group for the same proteins (first accession of a group; isoform stripped
# by DIA-NN's own inference here, kept for P04637-2 to exercise the normaliser)
DIANN = ["P02768", "P60709;P63261", "A0A075B6S2", "P04637", "Q11111", "Cont_P00761"]


def write_combined_protein(path):
    header = ["Protein", "Protein ID", "Entry Name", "Gene", "Protein Length", "Organism",
              "Protein Existence", "Description", "Protein Probability",
              "Top Peptide Probability", "Combined Total Peptides"]
    header += [f"{r} Spectral Count" for r in RUNS] + [f"{r} Intensity" for r in RUNS]
    header += [f"{r} MaxLFQ Intensity" for r in RUNS] + ["Indistinguishable Proteins"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for i, (prot, p_id, entry, gene) in enumerate(FRAGPIPE):
            lfq = [str(1000.0 * (i + 1) * (j + 1)) for j in range(len(RUNS))]
            w.writerow([prot, p_id, entry, gene, "500", "Homo sapiens", "1", "x", "1.0", "0.99",
                        "5"] + ["3"] * len(RUNS) + lfq + lfq + [""])


def write_report_tsv(path, groups):
    """A DE-contract report as TSV (compare_searches.py reads .tsv as well as parquet)."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Run", "Protein.Group", "PG.MaxLFQ", "Lib.PG.Q.Value"])
        for i, g in enumerate(groups):
            for j, r in enumerate(RUNS):
                w.writerow([r, g, 1000.0 * (i + 1) * (j + 1) * 1.1, 0.001])


def compare(out, *searches):
    """compare_searches.py on (label, report) pairs; its JSON."""
    args = [a for lab, path in searches for a in ("--search", f"{lab}:{path}")]
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "compare_searches.py"),
                        "--out", out, *args], capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise AssertionError(r.stderr)
    return json.loads(r.stdout)


class TheReading(unittest.TestCase):
    def test_accessions(self):
        cases = {
            "P12345;Q67890": "P12345",                 # DIA-NN group
            "P12345": "P12345",                        # FragPipe Protein ID
            "sp|P12345|ALBU_HUMAN": "P12345",          # FragPipe Protein, Sage
            "tr|A0A075B6S2|A0A075B6S2_HUMAN": "A0A075B6S2",   # 10-character accession
            "sp|P12345|A_HUMAN;sp|Q67890|B_HUMAN": "P12345",  # Sage group: the first member
            "P04637-2": "P04637",                      # isoform
            "sp|P04637-2|P53_HUMAN": "P04637",
            "contam_sp|P00167|CYB5_HUMAN": "contam_P00167",
            "sp|Cont_P00761|TRYP_PIG": "Cont_P00761",
            "Cont_P00761": "Cont_P00761",
            " P12345 ; Q1": "P12345",
            "": "", None: "",
        }
        for given, want in cases.items():
            with self.subTest(given=given):
                self.assertEqual(pid.normalize_protein_id(given), want)

    def test_fragpipe_protein_id(self):
        for prot, p_id, _, _ in FRAGPIPE:
            with self.subTest(prot=prot):
                got = pid.fragpipe_protein_id(prot, p_id)
                want = "contam_P00167" if prot.startswith("contam_") else p_id
                self.assertEqual(got, want)
        self.assertEqual(pid.fragpipe_protein_id("sp|P1|X_HUMAN", ""), "P1",
                         "no Protein ID: read the accession out of Protein")
        self.assertEqual(pid.fragpipe_protein_id("", "P1"), "P1")

    @unittest.skipUnless(shutil.which("Rscript"), "needs Rscript")
    def test_agrees_with_the_de_comparator_on_diann_and_fragpipe_ids(self):
        """compare_analyses.R's normalize_protein_id() compares the DE tables; the two must
        agree on what DIA-NN and FragPipe (Protein ID) write."""
        with open(os.path.join(SCRIPTS, "compare_analyses.R"), encoding="utf-8") as fh:
            src = fh.read()
        m = re.search(r"^normalize_protein_id <- function\(ids\) \{.*?^\}", src, re.S | re.M)
        self.assertIsNotNone(m)
        ids = ["P12345;Q67890", "P12345", "P04637-2", "A0A075B6S2", "Cont_P00761",
               "sp|P12345|ALBU_HUMAN", "sp|P04637-2|P53_HUMAN", "Q9Y6K9;Q9Y6K9-3"]
        # a file, not `Rscript -e`: -e mangles the function's backslashes on some platforms
        with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as fh:
            fh.write(m.group(0) + "\ncat(normalize_protein_id(c(" +
                     ", ".join(json.dumps(i) for i in ids) + ")), sep = '\\n')\n")
        try:
            r = subprocess.run(["Rscript", fh.name], capture_output=True, text=True, timeout=120)
        finally:
            os.unlink(fh.name)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stdout.split("\n")[:len(ids)],
                         [pid.normalize_protein_id(i) for i in ids])


@unittest.skipUnless(HAVE_ARROW, "needs pyarrow")
class FragPipeAgainstDiann(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.fp_out = os.path.join(self.d, "fragpipe")
        os.makedirs(self.fp_out)
        write_combined_protein(os.path.join(self.fp_out, "combined_protein.tsv"))
        self.diann = os.path.join(self.d, "diann_report.tsv")
        write_report_tsv(self.diann, DIANN)

    def tearDown(self):
        self._tmp.cleanup()

    def adapt(self):
        import pyarrow.parquet as pq
        import run_search
        report = run_search.adapt_fragpipe_dda(self.fp_out)
        return pq.read_table(report).to_pydict(), report

    def test_the_adapter_writes_accessions_and_keeps_the_names(self):
        t, _ = self.adapt()
        groups = list(dict.fromkeys(t["Protein.Group"]))
        self.assertEqual(groups, ["P02768", "P60709", "A0A075B6S2", "P04637-2", "Q99999",
                                  "contam_P00167", "Cont_P00761"])
        self.assertFalse(any("|" in g for g in groups))
        first = {g: (gn, nm) for g, gn, nm in zip(t["Protein.Group"], t["Genes"],
                                                  t["Protein.Names"])}
        self.assertEqual(first["P02768"], ("ALB", "ALBU_HUMAN"))
        self.assertEqual(first["contam_P00167"], ("CYB5A", "CYB5_HUMAN"))
        self.assertEqual(sorted(set(t["Run"])), list(RUNS))
        self.assertEqual(len(t["Run"]), len(FRAGPIPE) * len(RUNS))

    def test_compare_searches_finds_the_shared_proteins(self):
        _, report = self.adapt()
        j = compare(os.path.join(self.d, "cmp"), ("DIA-NN", self.diann), ("FragPipe", report))
        (pair,) = j["pairs"]
        # P02768, P60709, A0A075B6S2, P04637, Cont_P00761; contam_P00167 and Q99999 FragPipe-only,
        # Q11111 DIA-NN-only
        self.assertEqual(pair["shared"], 5, pair)
        self.assertEqual((pair["only_DIA-NN"], pair["only_FragPipe"]), (1, 2))
        self.assertEqual(pair["runs_compared"], 3)
        self.assertIsNotNone(pair["median_pearson_log2"])

    def test_a_report_adapted_before_the_fix_still_meets_diann(self):
        """Protein.Group = FragPipe's `Protein` (sp|...|...), as 2.8 and older wrote it."""
        old = os.path.join(self.d, "old_fragpipe.tsv")
        write_report_tsv(old, [p for p, _, _, _ in FRAGPIPE])
        j = compare(os.path.join(self.d, "cmp_old"), ("DIA-NN", self.diann), ("FragPipe", old))
        self.assertEqual(j["pairs"][0]["shared"], 5, j["pairs"][0])

    def test_two_groups_on_one_accession_are_counted(self):
        both = os.path.join(self.d, "isoforms.tsv")
        write_report_tsv(both, ["P04637", "P04637-2", "P02768"])
        j = compare(os.path.join(self.d, "cmp_iso"), ("A", both), ("B", self.diann))
        a = j["searches"][0]
        self.assertEqual((a["protein_groups"], a["rows_sharing_an_accession"]), (2, 3))
        self.assertEqual(a["accessions_shared_by_groups"], ["P04637"])
        self.assertEqual(j["searches"][1]["accessions_shared_by_groups"], [])
        with open(os.path.join(self.d, "cmp_iso", "SEARCH_COMPARISON.md")) as fh:
            md = fh.read()
        self.assertIn("3 rows named an accession another group", md)
        self.assertIn("kept (P04637).", md)


if __name__ == "__main__":
    unittest.main(verbosity=2)
