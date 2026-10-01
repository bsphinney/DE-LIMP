#!/usr/bin/env python3
"""The machine the DE ran on is on the record.

DPC-Quant's numbers are exact only on the same CPU family (PROT_0756 v2, zen4 vs zen2 HIVE
nodes: |dlogFC| <= 0.0022, |dt| <= 0.0064, no call changed). So:
  * run_de.R records `compute` in de_provenance.json -- CPU model, BLAS/LAPACK libraries,
    OPENBLAS_CORETYPE as found (never set) and, on a SLURM node, the node and its CPU-family
    feature (read with scontrol) -- and the same at the top of sessionInfo.txt;
  * provenance.py's REPRODUCE.md names that family, gives `sbatch --constraint=<family>` for an
    exact HIVE re-run, and warns against pinning OPENBLAS_CORETYPE; an older record without
    `compute` is said to be one.

A stub scontrol stands in for SLURM. Skips without R (arrow, dplyr, limma, jsonlite).
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
sys.path.insert(0, HERE)

from test_run_de_contaminants import SYNTH_R, r_has, rscript   # noqa: E402

NEEDS = ("limma", "arrow", "dplyr", "tidyr", "jsonlite")
STUB_SCONTROL = """#!/bin/sh
echo "NodeName=$3 Arch=x86_64 CoresPerSocket=56"
echo "   AvailableFeatures=cpu,ib,zen4"
echo "   ActiveFeatures=cpu,ib,zen4"
"""


@unittest.skipUnless(r_has(*NEEDS), "needs R with " + "/".join(NEEDS))
class RunDeRecordsCompute(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        bin_dir = os.path.join(cls.tmp, "bin")
        os.makedirs(bin_dir)
        with open(os.path.join(bin_dir, "scontrol"), "w") as fh:
            fh.write(STUB_SCONTROL)
        os.chmod(os.path.join(bin_dir, "scontrol"), 0o755)
        env = {k: v for k, v in os.environ.items() if k not in ("OPENBLAS_CORETYPE", "SLURMD_NODENAME")}
        slurm = dict(env, PATH=bin_dir + os.pathsep + env.get("PATH", ""), SLURMD_NODENAME="hive-dc-7-5-46")
        cls.runs = {}
        for key, e in (("plain", env), ("slurm", slurm)):
            out = os.path.join(cls.tmp, key)
            p = subprocess.run(["Rscript", os.path.join(SCRIPTS, "run_de.R"), "--input", "report.parquet",
                                "--metadata", "conditions.csv", "--method", "maxlfq", "--outdir", out],
                               capture_output=True, text=True, cwd=cls.tmp, env=e)
            cls.runs[key] = (p, out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def record(self, key):
        p, out = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            comp = json.load(fh)["compute"]
        with open(os.path.join(out, "sessionInfo.txt")) as fh:
            si = fh.read()
        return comp, si

    def test_the_machine_is_recorded(self):
        comp, si = self.record("plain")
        for k in ("cpu_model", "blas", "lapack", "lapack_version", "openblas_coretype", "note"):
            self.assertTrue(comp.get(k), (k, comp))
        self.assertEqual(comp["openblas_coretype"], "not set")     # recorded as found, never set
        self.assertIsNone(comp.get("cpu_family"))                   # no SLURM here
        self.assertNotIn("host", comp)                              # a laptop's name is not needed
        self.assertIn(f"CPU      {comp['cpu_model']}", si)
        self.assertIn(f"BLAS     {comp['blas']}", si)
        self.assertIn("OPENBLAS_CORETYPE not set", si)

    def test_a_slurm_node_names_its_cpu_family(self):
        comp, si = self.record("slurm")
        self.assertEqual(comp["slurm_node"], "hive-dc-7-5-46")
        self.assertEqual(comp["slurm_features"], ["cpu", "ib", "zen4"])
        self.assertEqual(comp["cpu_family"], "zen4")
        self.assertIn("(SLURM feature zen4, node hive-dc-7-5-46)", si)


class ReproduceMdNamesTheFamily(unittest.TestCase):
    def run_provenance(self, prov):
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        de = os.path.join(tmp, "de")
        os.makedirs(de)
        with open(os.path.join(de, "de_provenance.json"), "w") as fh:
            json.dump(prov, fh)
        out = os.path.join(tmp, "out")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"), "--outdir", out,
                            "--de-dir", de, "--timestamp", "2026-09-28T00:00:00Z"],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(os.path.join(out, "REPRODUCE.md"), encoding="utf-8") as fh:
            return fh.read()

    def test_family_constraint_and_the_coretype_warning(self):
        md = self.run_provenance({"method": "dpc", "compute": {
            "cpu_model": "AMD EPYC 9734 112-Core Processor", "slurm_node": "hive-dc-7-5-46",
            "cpu_family": "zen4", "blas": "libopenblasp-r0.3.34.so", "lapack": "libopenblasp-r0.3.34.so",
            "lapack_version": "3.12.0", "openblas_coretype": "not set"}})
        self.assertIn("## Exact numbers need the same kind of CPU", md)
        self.assertIn("`sbatch --constraint=zen4 ...`", md)
        self.assertIn("HIVE node `hive-dc-7-5-46` (CPU family `zen4`)", md)
        self.assertIn("|ΔlogFC| ≤ 0.0022 and |Δt| ≤ 0.0064", md)
        self.assertIn("Do not set\n`OPENBLAS_CORETYPE`", md)

    def test_without_a_family_it_says_how_to_find_one(self):
        md = self.run_provenance({"method": "dpc", "compute": {"cpu_model": "Apple M4", "blas": "x",
                                                               "lapack": "y"}})
        self.assertIn("`sbatch --constraint=<family>`", md)
        self.assertNotIn("--constraint=zen", md)

    def test_an_older_record_says_so(self):
        md = self.run_provenance({"method": "dpc"})
        self.assertIn("does not name the machine it ran on", md)


if __name__ == "__main__":
    unittest.main()
