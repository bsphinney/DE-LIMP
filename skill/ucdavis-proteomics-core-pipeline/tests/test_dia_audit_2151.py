#!/usr/bin/env python3
"""
DIA-route fixes from the 2026-11-17 audit of skill 2.15.0 (lane 1, and lane 2 #7). One test
class per finding; synthetic inputs only.

  L1-1  a Radiant library was reused whatever FASTA it was predicted from
  L1-2  the Radiant library was predicted with DIA-NN's defaults, not the skill's digest
  L1-3  the Radiant route's extraction widths are unvalidated: its record says experimental
  L1-5  the FragPipe-DIA adapter dropped the Global.* q-values and wrote a missing one as 0.0
  L1-7  the Radiant adapter wrote q = 0.0 when the output had no q-value column
  L1-11 a typed resolution overrode a timsTOF / Astral name
  L1-12 a DIA-NN cfg made without --fasta-meta had no --cont-quant-exclude
  L2-7  the de novo species call was not labelled exploratory
"""
import csv
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import estimate_params as ep  # noqa: E402

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_ARROW = True
except ImportError:                          # pragma: no cover
    HAVE_ARROW = False


def _exe(path, body):
    with open(path, "w") as fh:
        fh.write(body)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


# A DIA-NN stand-in for make_radiant_library.py: logs its argv, writes the predicted library for
# a prediction, and a parquet with Radiant's required columns for the conversion.
FAKE_DIANN = r"""#!/usr/bin/env python3
import os, sys
args = sys.argv[1:]
with open(os.environ["ARGV_LOG"], "a") as fh:
    fh.write(" ".join(args) + "\n")
outlib = args[args.index("--out-lib") + 1]
if "--fasta-search" in args:
    with open(outlib.rsplit(".", 1)[0] + ".predicted.speclib", "wb") as fh:
        fh.write(b"x" * 4000)
else:
    import pyarrow as pa, pyarrow.parquet as pq
    n = 40000   # the converted TSV must pass the reuse size floor (1 MB)
    pq.write_table(pa.table({
        "Precursor.Mz": [500.0 + i * 0.01 for i in range(n)], "Product.Mz": [300.0] * n,
        "Relative.Intensity": [1.0] * n, "Modified.Sequence": ["PEPTIDEK"] * n,
        "Precursor.Charge": [2] * n, "Fragment.Type": ["y"] * n, "Fragment.Charge": [1] * n,
        "Fragment.Series.Number": [3] * n}), outlib.rsplit(".", 1)[0] + ".parquet")
"""


@unittest.skipUnless(HAVE_ARROW, "needs pyarrow")
class RadiantLibrary(unittest.TestCase):
    """L1-1 / L1-2"""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.diann = _exe(os.path.join(self.d, "diann"), FAKE_DIANN)
        self.fasta = os.path.join(self.d, "db.fasta")
        with open(self.fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDEK\n")
        self.lib = os.path.join(self.d, "lib")
        self.log = os.path.join(self.d, "argv.txt")

    def build(self, *extra):
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_radiant_library.py"),
                            "--diann", self.diann, "--fasta", self.fasta, "--out-dir", self.lib,
                            *extra], capture_output=True, text=True, timeout=120,
                           env=dict(os.environ, ARGV_LOG=self.log))
        self.assertEqual(p.returncode, 0, p.stderr)
        return json.loads(p.stdout[p.stdout.index("{"):]), p

    def predictions(self):
        with open(self.log) as fh:
            return [ln.split() for ln in fh if "--fasta-search" in ln]

    def test_the_library_is_predicted_with_the_dia_nn_routes_digest(self):
        info, _ = self.build()
        (argv,) = self.predictions()
        line = " ".join(argv)
        self.assertIn(" ".join(ep.diann_library_flags()), line)
        for flag, want in (("--cut", ep.DIANN_DIGEST["cut"]), ("--min-pr-charge", "2"),
                           ("--max-pr-charge", "4"), ("--missed-cleavages", "1")):
            self.assertEqual(argv[argv.index(flag) + 1], want)
        self.assertIn("--met-excision", argv)
        self.assertIn("--unimod4", argv)
        self.assertEqual(info["library_flags"], ep.diann_library_flags() + ["--unimod4"])
        self.assertIn("DIA-NN's own default precursor m/z range", info["precursor_mz"])

    def test_the_same_inputs_reuse_it_and_a_changed_fasta_rebuilds_it(self):
        first, _ = self.build()
        self.assertNotIn("reused", first)
        with open(os.path.join(self.lib, "radiant_library.tsv.meta.json")) as fh:
            meta = json.load(fh)
        self.assertEqual(len(meta["fasta_sha256"]), 64)
        again, _ = self.build()
        self.assertTrue(again["reused"])
        self.assertEqual(len(self.predictions()), 1)
        with open(self.fasta, "a") as fh:            # a corrected database, same path
            fh.write(">sp|P2|Y\nANOTHERK\n")
        rebuilt, p = self.build()
        self.assertNotIn("reused", rebuilt)
        self.assertIn("fasta_sha256", p.stderr)
        self.assertEqual(len(self.predictions()), 2)

    def test_the_build_flags_are_the_cfgs(self):
        """One definition: the cfg's own digest/charge/fragment lines are the library flags."""
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                                "--engine", "diann", "--acquisition", "DIA", "--instrument",
                                "timsTOF HT", "--precursor-mz-range", "300", "1200", "--out", cfg],
                               capture_output=True, text=True, timeout=60)
            self.assertEqual(p.returncode, 0, p.stderr)
            lines = open(cfg).read().split()
        flags = ep.diann_library_flags()
        for i, tok in enumerate(flags):
            if tok.startswith("--"):
                self.assertIn(tok, lines)
                if i + 1 < len(flags) and not flags[i + 1].startswith("--"):
                    self.assertEqual(lines[lines.index(tok) + 1], flags[i + 1], tok)


class RadiantPreset(unittest.TestCase):
    """L1-3"""

    def test_the_preset_record_says_experimental(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "r.radiantConfig")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_presets.py"),
                                "--engine", "radiant", "--instrument", "Orbitrap Astral",
                                "--out", out], capture_output=True, text=True, timeout=60)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out + ".provenance.json") as fh:
                prov = json.load(fh)
        self.assertEqual(prov["route_status"]["status"], "experimental")
        self.assertIn("not been benchmarked", prov["route_status"]["reasons"][0])
        self.assertTrue(any("UNVALIDATED" in v for v in prov["keys_changed"].values()))


@unittest.skipUnless(HAVE_ARROW, "needs pyarrow")
class Adapters(unittest.TestCase):
    """L1-5 / L1-7"""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        import run_search
        self.rs = run_search

    def radiant(self, cols):
        d = os.path.join(self.d, "radiant_results", "fulcrum-results")
        os.makedirs(d)
        pq.write_table(pa.table(cols), os.path.join(d, "part-00000.parquet"))

    def test_radiant_without_a_q_value_column_is_refused(self):
        self.radiant({"Run": ["a.mzML", "b.mzML"], "Protein.Group": ["P1", "P1"],
                      "PG.Normalised": [1.0, 2.0]})
        with self.assertRaises(SystemExit) as cm:
            self.rs.adapt_radiant(self.d)
        self.assertIn("no q-value column", str(cm.exception))
        self.assertFalse(os.path.exists(os.path.join(self.d, "report.parquet")))
        with self.assertRaises(ValueError):
            self.rs.radiant_quantity_declaration(self.d, "PG.Normalised", None)

    def test_radiant_rows_without_a_q_value_are_left_out_never_given_zero(self):
        self.radiant({"Run": ["a.mzML", "a.mzML", "b.mzML", "b.mzML"],
                      "Protein.Group": ["P1", "P1", "P2", "P1"],
                      "PG.Normalised": [5.0, 3.0, 2.0, 4.0],
                      "Global.PG.Q.Value": [None, 0.02, None, 0.001]})
        rep = self.rs.adapt_radiant(self.d)
        t = pq.read_table(rep).to_pydict()
        rows = {(r, p): (v, q) for r, p, v, q in zip(t["Run"], t["Protein.Group"],
                                                     t["PG.MaxLFQ"], t["Q.Value"])}
        # P1 in a: highest quantity, best q among the rows that have one
        self.assertEqual(rows[("a", "P1")], (5.0, 0.02))
        self.assertEqual(rows[("b", "P1")], (4.0, 0.001))
        self.assertNotIn(("b", "P2"), rows, "a protein x run with no q-value was given one")

    def test_a_precursor_q_value_is_not_called_protein_level(self):
        dec = self.rs.radiant_quantity_declaration(self.d, "PG.Quantity", "Q.Value")
        self.assertIn("PRECURSOR q-value", dec["identification"])
        self.assertIn("not a protein-level FDR", dec["identification"])
        dec = self.rs.radiant_quantity_declaration(self.d, "PG.Quantity", "Global.PG.Q.Value")
        self.assertNotIn("PRECURSOR", dec["identification"])

    def fragpipe_tsv(self, cols, rows):
        os.makedirs(os.path.join(self.d, "dia-quant-output"))
        with open(os.path.join(self.d, "dia-quant-output", "report.tsv"), "w") as fh:
            fh.write("\t".join(cols) + "\n" + "".join("\t".join(r) + "\n" for r in rows))

    def test_fragpipe_dia_carries_the_global_q_values(self):
        cols = ["Run", "Protein.Group", "PG.MaxLFQ", "Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value",
                "Global.Q.Value", "Global.PG.Q.Value", "PG.Q.Value"]
        self.fragpipe_tsv(cols, [["S1", "P1", "10", "0.001", "0.002", "0.003", "0.004",
                                  "0.05", "0.006"],
                                 ["S2", "P1", "12", "0.001", "0.002", "0.003", "", "0.2", ""]])
        t = pq.read_table(self.rs.adapt_fragpipe_dia(self.d)).to_pydict()
        self.assertEqual(t["Global.PG.Q.Value"], [0.05, 0.2])
        self.assertEqual(t["Global.Q.Value"][0], 0.004)
        self.assertNotEqual(t["Global.Q.Value"][1], t["Global.Q.Value"][1], "empty -> NaN")
        self.assertEqual(t["PG.Q.Value"][0], 0.006)

    def test_fragpipe_dia_without_a_required_q_column_is_refused(self):
        self.fragpipe_tsv(["Run", "Protein.Group", "PG.MaxLFQ", "Q.Value"],
                          [["S1", "P1", "10", "0.001"]])
        with self.assertRaises(SystemExit) as cm:
            self.rs.adapt_fragpipe_dia(self.d)
        self.assertIn("Lib.Q.Value", str(cm.exception))
        self.assertFalse(os.path.exists(os.path.join(self.d, "report.parquet")))


class InstrumentName(unittest.TestCase):
    """L1-11"""

    def test_an_astral_keeps_the_readmes_ms2_value_whatever_the_resolution(self):
        cls, ms1, ms2, _, _ = ep.classify_instrument("Orbitrap Astral", 240000, 80000)
        self.assertEqual((cls, ms1, ms2), ("orbitrap_astral", 4, 10))
        cls, ms1, ms2, label, _ = ep.classify_instrument("Orbitrap Astral", 180000, 80000)
        self.assertEqual((cls, ms1, ms2), ("orbitrap_astral", 5.0, 10))
        self.assertIn("180,000", label)

    def test_a_timstof_with_resolutions_is_still_a_timstof(self):
        cls, ms1, ms2, _, src = ep.classify_instrument("timsTOF HT", 60000, 30000)
        self.assertEqual((cls, ms1, ms2), ("timstof", 15, 15))
        self.assertIn("not used", src)
        self.assertEqual(ep.classify_instrument("ZenoTOF 7600", 30000, 30000)[:3],
                         ("sciex_tof", 20, 20))

    def test_an_orbitrap_resolution_still_wins_for_an_orbitrap(self):
        self.assertEqual(ep.classify_instrument("Orbitrap Exploris 480", 120000, 30000)[:3],
                         ("orbitrap_measured", 7, 15))


class FastaMetaDefault(unittest.TestCase):
    """L1-12"""

    def run_in(self, d, *extra):
        cfg = os.path.join(d, "wf", "p.cfg")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "diann", "--acquisition", "DIA", "--instrument",
                            "timsTOF HT", "--precursor-mz-range", "300", "1200", "--out", cfg,
                            *extra], capture_output=True, text=True, timeout=60, cwd=d)
        self.assertEqual(p.returncode, 0, p.stderr)
        return open(cfg).read().splitlines(), p.stderr

    def test_the_sessions_sidecar_is_used_when_none_is_given(self):
        with tempfile.TemporaryDirectory() as d:
            lines, err = self.run_in(d)
            self.assertFalse([ln for ln in lines if ln.startswith("--cont-quant-exclude")])
            self.assertIn("WARNING: no FASTA sidecar", err)
            with open(os.path.join(d, "search.fasta.meta.json"), "w") as fh:
                json.dump({"diann_cont_quant_exclude": "Cont_"}, fh)
            lines, err = self.run_in(d)
            self.assertIn("--cont-quant-exclude Cont_", lines)
            self.assertIn("--fasta-meta not given: using", err)


class DenovoExploratory(unittest.TestCase):
    """L2-7"""

    def test_the_species_call_says_exploratory_and_closest_relative(self):
        with tempfile.TemporaryDirectory() as d:
            # a tiny taxonomy: 1 root > 33208 Metazoa > 9605 Homo (genus) > 9606 H. sapiens
            nodes = os.path.join(d, "nodes.dmp")
            names = os.path.join(d, "names.dmp")
            with open(nodes, "w") as fh:
                for t, p_, r in ((1, 1, "no rank"), (33208, 1, "kingdom"), (9605, 33208, "genus"),
                                 (9606, 9605, "species")):
                    fh.write(f"{t}\t|\t{p_}\t|\t{r}\t|\t\t|\n")
            with open(names, "w") as fh:
                for t, n in ((9606, "Homo sapiens"), (9605, "Homo"), (33208, "Metazoa")):
                    fh.write(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n")
            hits = os.path.join(d, "hits.tsv")
            with open(hits, "w") as fh:
                for pep in ("PEPTIDEA", "PEPTIDEB"):
                    for _ in range(3):
                        fh.write("\t".join([pep, "s", "100", "8", "0", "0", "1", "8", "1", "8",
                                            "1e-3", "30", "9606"]) + "\n")
            prefix = os.path.join(d, "out")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "denovo_lca.py"), hits,
                                prefix, "--nodes", nodes, "--names", names],
                               capture_output=True, text=True, timeout=60)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertTrue(p.stdout.startswith("EXPLORATORY:"), p.stdout)
            self.assertIn("closest sequenced relatives", p.stdout)
            self.assertNotIn("top HOST species", p.stdout)
            with open(prefix + "_species_call.json") as fh:
                rec = json.load(fh)
            self.assertEqual(rec["status"], "exploratory")
            self.assertEqual(rec["closest_sequenced_relatives_host"],
                             [{"taxon": "Homo sapiens (species)", "peptides": 2}])
            with open(prefix + "_peptide_lca.tsv") as fh:
                self.assertEqual(len(list(csv.reader(fh, delimiter="\t"))), 3)

    def test_the_fdr_record_never_calls_its_cut_a_sequence_fdr(self):
        import denovo_homology_fdr as dh
        self.assertIn("NOT a false discovery rate on de novo sequence correctness",
                      dh.FDR_MEANING)
        with open(os.path.join(SCRIPTS, "denovo_homology_fdr.py")) as fh:
            src = fh.read()
        self.assertIn('"fdr_meaning": FDR_MEANING', src)
        self.assertIn('"status": EXPLORATORY["status"]', src)
        with open(os.path.join(os.path.dirname(SCRIPTS), "references", "denovo.md")) as fh:
            ref = fh.read()
        self.assertIn("Exploratory, not a species identification", ref)
        self.assertNotIn("# 3. FDR", ref)


if __name__ == "__main__":
    unittest.main()
