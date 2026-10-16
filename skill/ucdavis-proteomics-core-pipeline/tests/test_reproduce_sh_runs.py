#!/usr/bin/env python3
"""
reproduce.sh (provenance.py) is RUN here, not only read: steps 5 and 6 execute in bash with stub
run_search.py / Rscript that record the arguments they were given.

Why: skill 2.13.0's SILAC note was inserted after the trailing backslash of step 5's
`run_search.py --tools ... \\` line. A comment there ends the continuation, so bash ran
run_search.py with no --bundle/--params/--fasta/--out and then `--bundle` as a command of its own
("--bundle: command not found", exit 127). `bash -n` and substring checks both passed.

Pinned:
  * a SILAC bundle: step 5 searches with the recorded cfg and every argument; step 6 replays
    run_silac.R with the bundle's design and cfg;
  * a maxlfq + raw bundle (its no-norm note sits on step 5 too): every argument still arrives;
  * no generated reproduce.sh has a line ending in a backslash followed by a comment or a blank;
  * the DE's design is copied into inputs/ without replacing a different file of the same name,
    and the same file reached by a symlink is not copied twice.
Stdlib only.
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

STUB = """#!/usr/bin/env python3
import json, os, sys
with open(os.environ["ARGV_LOG"], "a") as fh:
    fh.write(json.dumps({"prog": %r, "argv": sys.argv[1:]}) + "\\n")
"""


def provenance(out, de, *more):
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"), "--outdir", out,
                        "--de-dir", de, *more], capture_output=True, text=True, timeout=300)
    assert p.returncode == 0, p.stdout[-2000:] + p.stderr[-2000:]
    with open(os.path.join(out, "reproduce.sh")) as fh:
        return fh.read()


def continuation_breaks(sh):
    """Lines ending in a backslash whose next line is a comment or blank: bash ends the command
    there and runs the rest as another command."""
    lines = sh.splitlines()
    return [(i + 1, lines[i], lines[i + 1]) for i in range(len(lines) - 1)
            if lines[i].rstrip().endswith("\\") and
            (not lines[i + 1].strip() or lines[i + 1].lstrip().startswith("#"))]


def run_steps_5_and_6(d, out, sh):
    """Execute reproduce.sh from step 5 to the end with stubs; -> the calls, in order."""
    skill = os.path.join(d, "skill")
    os.makedirs(os.path.join(skill, "scripts"), exist_ok=True)
    bindir = os.path.join(d, "bin")
    os.makedirs(bindir, exist_ok=True)
    for path, prog in ((os.path.join(skill, "scripts", "run_search.py"), "run_search.py"),
                       (os.path.join(bindir, "Rscript"), "Rscript")):
        with open(path, "w") as fh:
            fh.write(STUB % prog)
        os.chmod(path, 0o755)
    body = sh[sh.index("\n# 5. "):sh.index('\necho "Done.')]
    log = os.path.join(d, "argv.jsonl")
    env = dict(os.environ, SKILL=skill, ARGV_LOG=log,
               PATH=bindir + os.pathsep + os.environ["PATH"])
    p = subprocess.run(["bash", "-c", "set -euo pipefail\n" + body], cwd=out, env=env,
                       capture_output=True, text=True, timeout=120)
    assert p.returncode == 0, p.stderr
    with open(log) as fh:
        return [json.loads(ln) for ln in fh], p.stderr


def arg(argv, flag):
    return argv[argv.index(flag) + 1] if flag in argv else None


class SilacBundle(unittest.TestCase):
    def make(self, d, conditions=None, design_name="silac_design.csv", design_text=None):
        de = os.path.join(d, "tables")
        os.makedirs(de)
        design = os.path.join(d, "analysis", design_name)
        os.makedirs(os.path.dirname(design))
        with open(design, "w") as fh:
            fh.write(design_text or "File.Name,Replicate,Channel,Group\nr1,R1,L,WT\nr1,R1,H,KO\n")
        cfg = os.path.join(d, "params.cfg")
        with open(cfg, "w") as fh:
            fh.write("--fixed-mod SILAC,0.0,KR,label\n--channels SILAC,L,KR,0:0;SILAC,H,KR,8:10\n")
        with open(cfg + ".rationale.json", "w") as fh:
            json.dump({"labels": {"scheme": "silac2", "channels": []}}, fh)
        with open(os.path.join(de, "de_provenance.json"), "w") as fh:
            json.dump({"pipeline_id": "silac_ratio", "method": "silac", "metadata": design,
                       "replay": {"script": "run_silac.R", "args": [
                           "--input", "/elsewhere/report.parquet", "--design", design,
                           "--outdir", "/elsewhere/tables", "--params", cfg,
                           "--comparisons", "KO-WT", "--run-level-only", "--adjp", "0.05"]}}, fh)
        out = os.path.join(d, "repro")
        more = ["--params", cfg, "--engine", "diann", "--acquisition", "DIA",
                "--instrument", "timsTOF HT"]
        if conditions:
            more += ["--conditions", conditions]
        return out, provenance(out, de, *more), design

    def test_steps_5_and_6_run_with_every_argument(self):
        with tempfile.TemporaryDirectory() as d:
            out, sh, _ = self.make(d)
            self.assertEqual(continuation_breaks(sh), [])
            calls, _ = run_steps_5_and_6(d, out, sh)
        self.assertEqual([c["prog"] for c in calls], ["run_search.py", "Rscript"])
        s = calls[0]["argv"]
        self.assertEqual(arg(s, "--bundle"), "./wf/workflow.manifest.json")
        self.assertEqual(arg(s, "--params"), "./inputs/params.cfg", "the recorded labels")
        self.assertEqual(arg(s, "--fasta"), "./search.fasta")
        self.assertEqual(arg(s, "--out"), "./search_out")
        self.assertIn("--files", s)
        r = calls[1]["argv"]
        self.assertTrue(r[0].endswith("/scripts/run_silac.R"), r)
        self.assertEqual(arg(r, "--input"), "./search_out/report.parquet")
        self.assertEqual(arg(r, "--design"), "inputs/silac_design.csv")
        self.assertEqual(arg(r, "--params"), "inputs/params.cfg")
        self.assertIn("--run-level-only", r)

    def test_a_different_file_of_the_same_name_is_never_replaced(self):
        with tempfile.TemporaryDirectory() as d:
            cond = os.path.join(d, "silac_design.csv")
            with open(cond, "w") as fh:
                fh.write("File.Name,Group\nr1,WT\n")
            out, sh, design = self.make(d, conditions=cond)
            with open(os.path.join(out, "inputs", "silac_design.csv")) as fh:
                self.assertEqual(fh.read(), "File.Name,Group\nr1,WT\n", "the conditions stay")
            with open(os.path.join(out, "inputs", "silac_design.de_design.csv")) as fh, \
                    open(design) as src:
                self.assertEqual(fh.read(), src.read())
            with open(os.path.join(out, "MANIFEST.txt")) as fh:
                self.assertIn("both kept", fh.read())
            calls, _ = run_steps_5_and_6(d, out, sh)
        self.assertEqual(arg(calls[1]["argv"], "--design"), "inputs/silac_design.de_design.csv")

    def test_the_same_file_by_another_path_is_not_copied_twice(self):
        with tempfile.TemporaryDirectory() as d:
            real = os.path.join(d, "real_design.csv")
            with open(real, "w") as fh:
                fh.write("File.Name,Replicate,Channel,Group\nr1,R1,L,WT\nr1,R1,H,KO\n")
            out, sh, design = self.make(d, design_name="link.csv")
            shutil.rmtree(out)                      # a fresh bundle: the symlinked design only
            os.unlink(design)
            os.symlink(real, design)
            de = os.path.join(d, "tables")
            sh = provenance(out, de, "--conditions", real, "--params",
                            os.path.join(d, "params.cfg"))
            self.assertEqual(sorted(f for f in os.listdir(os.path.join(out, "inputs"))
                                    if f.endswith(".csv")), ["real_design.csv"])
            calls, _ = run_steps_5_and_6(d, out, sh)
        self.assertEqual(arg(calls[1]["argv"], "--design"), "inputs/real_design.csv")


class LabelFreeBundle(unittest.TestCase):
    def test_a_maxlfq_raw_bundle_still_passes_every_search_argument(self):
        with tempfile.TemporaryDirectory() as d:
            de = os.path.join(d, "tables")
            os.makedirs(de)
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump({"method": "maxlfq", "quantities": "raw"}, fh)
            out = os.path.join(d, "repro")
            sh = provenance(out, de, "--de-method", "maxlfq")
            self.assertEqual(continuation_breaks(sh), [])
            calls, err = run_steps_5_and_6(d, out, sh)
        self.assertIn("re-quantify with --no-norm", err)
        s = calls[0]["argv"]
        for flag in ("--bundle", "--params", "--fasta", "--out", "--files"):
            self.assertIn(flag, s)
        self.assertTrue(calls[1]["argv"][0].endswith("/scripts/run_de.R"))
        self.assertEqual(arg(calls[1]["argv"], "--quantities"), "raw")


if __name__ == "__main__":
    unittest.main()
