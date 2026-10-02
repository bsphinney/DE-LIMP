#!/usr/bin/env python3
"""
The DIA-NN 5-step chain in a folder whose name has a space.

Core service folders hold spaces (2,557 directories within three levels of Data/lab/service).
The chain quoted the raw files but spliced --out, --fasta, --lib and --temp unquoted: step 2's
`sed -n ... <out>/file_list.txt` split on the space ("no file for task 0") and step 1 handed
DIA-NN half a FASTA path (HIVE-ops review of 2.10). Here every generated job -- steps 1, 1b, 2,
3, 4 and 5 -- runs in order against a fake DIA-NN that refuses any path it cannot open, so a path
split anywhere fails the test.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)
from job_env import job_env  # noqa: E402

FAKE = r"""#!/usr/bin/env python3
import json, os, sys
a = sys.argv[1:]
def val(f):
    return a[a.index(f) + 1] if f in a else None
with open(os.environ["FAKE_ARGV"], "a") as fh:
    fh.write(json.dumps(a) + "\n")
for p in (val("--fasta"), val("--lib")) + tuple(val("--f") for _ in [0] if "--f" in a):
    if p and not os.path.exists(p):
        print("ERROR: cannot open", repr(p)); sys.exit(1)
for f in ("--temp", "--out", "--out-lib"):
    d = os.path.dirname(val(f) or "") if f != "--temp" else val(f)
    if d and not os.path.isdir(d):
        print("ERROR: no folder", repr(d), "for", f); sys.exit(1)
if "--predictor" in a:
    lib = val("--out-lib")
    open(lib[:-len(".speclib")] + ".predicted.speclib", "w").write("lib")
    sys.exit(0)
print("[0:04] Scan window radius set to 7"); sys.stdout.flush()
if "--f" in a and val("--temp"):
    f = val("--f")
    open(os.path.join(val("--temp"), os.path.splitext(os.path.basename(f))[0] + ".quant"),
         "w").write("q")
if val("--out-lib"):
    open(val("--out-lib"), "w").write("lib")
if val("--out") and "--use-quant" in a:
    open(val("--out"), "w").write("report")
"""


class ChainInAFolderWithASpace(unittest.TestCase):
    def test_every_step_runs_with_its_paths_whole(self):
        with tempfile.TemporaryDirectory() as d:
            base = os.path.join(d, "Lab Name 09182015", "Service run")
            raws = []
            for i in range(2):
                raws.append(os.path.join(base, "raw data", f"run {i}.mzML"))
                os.makedirs(os.path.dirname(raws[-1]), exist_ok=True)
                with open(raws[-1], "wb") as fh:
                    fh.truncate((1 + i) * 1024 * 1024)
            fasta = os.path.join(base, "db v1.fasta")
            with open(fasta, "w") as fh:
                fh.write(">sp|P1|X\nPEPTIDEK\n")
            cfg = os.path.join(base, "params v1.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n")   # step 1b probes
            diann = os.path.join(d, "diann-linux")
            with open(diann, "w") as fh:
                fh.write(FAKE)
            os.chmod(diann, 0o755)
            out = os.path.join(base, "output", "search out")
            argv_log = os.path.join(d, "argv.jsonl")
            env = job_env(d, PATH=os.environ["PATH"], FAKE_ARGV=argv_log)
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", diann, "--raw", *raws, "--fasta", fasta, "--out", out,
                                "--cfg", cfg, "--partition", "low", "--account", "publicgrp"],
                               capture_output=True, text=True, env=env, timeout=120)
            self.assertEqual(g.returncode, 0, g.stderr)
            steps = [("step1_libpred", None), ("step1b_window", None),
                     ("step2_firstpass", "0"), ("step2_firstpass", "1"), ("step3_assembly", None),
                     ("step4_finalpass", "0"), ("step4_finalpass", "1"), ("step5_report", None)]
            for name, task in steps:
                e = dict(env, **({"SLURM_ARRAY_TASK_ID": task} if task else {}))
                # job_env: not a search job's side effects -- every hook is switched off in env
                p = subprocess.run(["bash", os.path.join(out, f"{name}.sbatch")], cwd=out,
                                   capture_output=True, text=True, env=e, timeout=240)
                self.assertEqual(p.returncode, 0, f"{name} {task}: {p.stdout}\n{p.stderr}")
            self.assertTrue(os.path.isfile(os.path.join(out, "report.parquet")))
            with open(argv_log) as fh:
                calls = [json.loads(ln) for ln in fh]
            self.assertTrue(calls)
            for c in calls:
                if "--fasta" in c:
                    self.assertEqual(c[c.index("--fasta") + 1], fasta)
            self.assertEqual(open(os.path.join(out, "window.txt")).read().strip(), "7")


if __name__ == "__main__":
    unittest.main()
