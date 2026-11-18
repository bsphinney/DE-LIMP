#!/usr/bin/env python3
"""
The 5-step chain's step scripts, run with bash as a compute node runs them, against a stand-in
DIA-NN. Audit 2026-11-17 (skill 2.15.0):

  lane 5 #1   a step judged success only by a non-empty file, never by DIA-NN's exit status: a
              stand-in that wrote 7 bytes of .quant and exited 139 left step 2 COMPLETED.
  lane 5 #10  step 3 copied quant_step2 to quant_step2_orig, which nothing read, and a second run
              nested the copy inside itself.
  lane 5 #16  the logs were relative (`#SBATCH -o s2_firstpass_%A_%a.log`), so a step
              resubmitted by path from $HOME logged to $HOME, where the watcher never looks.
"""
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
from job_env import job_env  # noqa: E402  (env for running job scripts)

# Writes what DIA-NN would (a .quant under --temp for one file, --out-lib / --out otherwise),
# then exits with $FAKE_RC.
FAKE_DIANN = r"""#!/bin/bash
f=""; tmp=""; out=""; outlib=""
while [ $# -gt 0 ]; do
  case "$1" in --f) f="$2"; shift;; --temp) tmp="$2"; shift;; --out) out="$2"; shift;;
    --out-lib) outlib="$2"; shift;; esac; shift
done
if [ -n "$f" ] && [ -n "$tmp" ]; then b="${f##*/}"; printf 'partial' > "$tmp/${b%.*}.quant"; fi
if [ -n "$outlib" ]; then printf 'partial-lib' > "$outlib"; fi
if [ -n "$out" ]; then mkdir -p "$(dirname "$out")"; printf 'partial-report' > "$out"; fi
echo "fake diann: exiting ${FAKE_RC:-0}"
exit "${FAKE_RC:-0}"
"""


def _exe(path, body):
    with open(path, "w") as fh:
        fh.write(body)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def _read(path):
    with open(path) as fh:
        return fh.read()


class ChainEngineStatus(unittest.TestCase):

    def _chain(self, d, n=3):
        raws = []
        for i in range(n):
            p = os.path.join(d, f"run{i}.mzML")
            with open(p, "w") as fh:
                fh.write("x")
            raws.append(p)
        fasta = os.path.join(d, "db.fasta")
        with open(fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDER\n")
        cfg = os.path.join(d, "p.cfg")
        with open(cfg, "w") as fh:
            # mass accuracy and window pinned: no step 1b
            fh.write("--qvalue 0.01\n--mass-acc 20\n--mass-acc-ms1 7\n--window 8\n")
        diann = _exe(os.path.join(d, "diann"), FAKE_DIANN)
        out = os.path.join(d, "out")
        g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                            "--diann", diann, "--raw", *raws, "--fasta", fasta, "--out", out,
                            "--cfg", cfg, "--threads-per-file", "4"],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(g.returncode, 0, g.stderr)
        return out, json.loads(g.stdout)

    def _run(self, d, out, name, rc, task=0):
        env = job_env(d, FAKE_RC=str(rc), SLURM_ARRAY_TASK_ID=str(task))
        return subprocess.run(["bash", os.path.join(out, name)], cwd=out, capture_output=True,
                              text=True, env=env, timeout=120)

    def test_a_crash_that_leaves_a_partial_quant_fails_step_2_and_the_file_goes(self):
        with tempfile.TemporaryDirectory() as d:
            out, _ = self._chain(d)
            quant = os.path.join(out, "quant_step2", "run0.quant")
            p = self._run(d, out, "step2_firstpass.sbatch", 139)
            log = p.stdout + p.stderr
            self.assertNotEqual(p.returncode, 0, log)
            self.assertIn("FAILED: DIA-NN exited with status 139", log)
            self.assertFalse(os.path.exists(quant), "the partial .quant was left for step 4")
            # the same engine exiting 0 passes, and its .quant stays
            p = self._run(d, out, "step2_firstpass.sbatch", 0)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertTrue(os.path.exists(quant))

    def test_every_dia_nn_step_reads_the_exit_status(self):
        with tempfile.TemporaryDirectory() as d:
            out, _ = self._chain(d)
            names = ["step1_libpred.sbatch", "step2_firstpass.sbatch", "step3_assembly.sbatch",
                     "step4_finalpass.sbatch", "step5_report.sbatch"]
            for name in names:
                body = _read(os.path.join(out, name))
                self.assertIn("|| DIANN_RC=$?", body, name)
                # the status is read before the artefact is
                self.assertLess(body.index('if [ "$DIANN_RC" -ne 0 ]'),
                                body.index("did not write"), name)
            # step 1 and step 3, run: a non-zero exit fails them and removes the library
            p = self._run(d, out, "step1_libpred.sbatch", 134)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("FAILED: DIA-NN exited with status 134", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(out, "step1.predicted.speclib")))

    def test_step_3_makes_no_quant_step2_orig_copy(self):
        with tempfile.TemporaryDirectory() as d:
            out, _ = self._chain(d)
            self.assertNotIn("quant_step2_orig", _read(os.path.join(out, "step3_assembly.sbatch")))

    def test_every_step_logs_into_the_search_folder_by_absolute_path(self):
        with tempfile.TemporaryDirectory() as d:
            out, _ = self._chain(d)
            for name in ("step1_libpred.sbatch", "step2_firstpass.sbatch",
                         "step3_assembly.sbatch", "step4_finalpass.sbatch",
                         "step5_report.sbatch"):
                body = _read(os.path.join(out, name))
                logs = re.findall(r'^#SBATCH -[oe] "(.+)"$', body, re.M)
                self.assertEqual(len(logs), 2, name)
                for log in logs:
                    self.assertTrue(os.path.isabs(log), (name, log))
                    self.assertEqual(os.path.dirname(log), os.path.abspath(out), (name, log))
            self.assertIn(f'#SBATCH -o "{os.path.abspath(out)}/s2_firstpass_%A_%a.log"',
                          _read(os.path.join(out, "step2_firstpass.sbatch")))


if __name__ == "__main__":
    unittest.main()
