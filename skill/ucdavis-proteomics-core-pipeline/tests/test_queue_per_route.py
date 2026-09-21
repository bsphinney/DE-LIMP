#!/usr/bin/env python3
"""
Which queue each job-array route actually writes into its #SBATCH headers -- the behaviour
references/environment.md ("SLURM submission") describes. These pin existing behaviour so
that the doc and the code cannot drift apart again without a test noticing.

The trap they pin: slurm_queue() has a rule for preemption-safe ARRAY steps (move to
publicgrp/low when fewer than 2 x need CPUs are free on `high`, or when `squeue` cannot run).
radiant_parallel.py asks slurm_queue() separately for its array, so the rule applies there.
diann_parallel.py does not: main() fills in --partition AND --account from one
slurm_queue() call (no peak, not preemption-safe) before it asks per step, and slurm_queue()
returns an explicit partition+account untouched -- so the per-step call for steps 2 and 4
never reaches the array rule, and all five steps share one queue. An earlier version of the
doc said steps 2 and 4 move to `low` on their own; the generated scripts said `high`.

If diann_parallel.py is changed to give its array steps their own queue, update
references/environment.md together with these tests.

Everything is stubbed: sacctmgr (a facility member who also has publicgrp/low), sinfo (5000
idle CPUs on low) and, per test, squeue. PATH holds only the stub directory.
"""
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")

SACCTMGR = "#!/bin/sh\necho 'genome-center-grp|high|'\necho 'publicgrp|low|publicgrp-low-qos'\n"
SINFO = "#!/bin/sh\necho '0/5000/0/5000'\n"


def _exe(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, 0o755)


def queue_of(path):
    """(partition, account, qos, requeue) from a generated script's #SBATCH header."""
    txt = open(path).read()

    def get(opt):
        m = re.search(rf"^#SBATCH --{opt}=(\S+)$", txt, re.M)
        return m.group(1) if m else None
    return (get("partition"), get("account"), get("qos"),
            bool(re.search(r"^#SBATCH --requeue$", txt, re.M)))


class QueueHarness(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.bin = os.path.join(self.d, "bin")
        os.makedirs(self.bin)
        _exe(os.path.join(self.bin, "sacctmgr"), SACCTMGR)
        _exe(os.path.join(self.bin, "sinfo"), SINFO)

    def tearDown(self):
        self._tmp.cleanup()

    def squeue_reports(self, cpus_in_use):
        """None: no squeue at all (slurm_queue() then cannot read your usage)."""
        if cpus_in_use is not None:
            _exe(os.path.join(self.bin, "squeue"), f"#!/bin/sh\necho {cpus_in_use}\n")

    def run_script(self, argv):
        env = {"PATH": self.bin, "HOME": self.d, "USER": "someone"}
        r = subprocess.run([sys.executable, *argv], cwd=self.d, capture_output=True,
                           text=True, env=env, timeout=120)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return r


class DiannChainQueueTests(QueueHarness):
    STEPS = ("step1_libpred", "step2_firstpass", "step3_assembly", "step4_finalpass",
             "step5_report")

    def chain(self):
        raws = []
        for i in range(6):
            raws.append(os.path.join(self.d, f"f{i}.d"))
            os.makedirs(raws[-1])
        cfg = os.path.join(self.d, "diann.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n--window 7\n")
        fasta = os.path.join(self.d, "db.fasta")
        with open(fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDEK\n")
        out = os.path.join(self.d, "out")
        self.run_script([os.path.join(SCRIPTS, "diann_parallel.py"), "--diann",
                         "/opt/diann-2.6.1/diann-linux", "--raw", *raws, "--fasta", fasta,
                         "--out", out, "--cfg", cfg, "--no-probe-window"])
        return {s: queue_of(os.path.join(out, f"{s}.sbatch")) for s in self.STEPS}

    def test_array_steps_do_not_move_to_low_when_squeue_cannot_run(self):
        self.squeue_reports(None)
        q = self.chain()
        for step, (part, acct, _qos, _rq) in q.items():
            self.assertEqual((part, acct), ("high", "genome-center-grp"), step)

    def test_array_steps_do_not_move_to_low_with_fewer_than_2x16_free(self):
        self.squeue_reports(40)                     # 24 of 64 free: < 2 x 16, >= 16
        q = self.chain()
        for step, (part, acct, _qos, _rq) in q.items():
            self.assertEqual((part, acct), ("high", "genome-center-grp"), step)

    def test_the_whole_chain_moves_to_low_when_fewer_than_16_are_free(self):
        self.squeue_reports(50)                     # 14 of 64 free
        q = self.chain()
        for step, (part, acct, qos, requeue) in q.items():
            self.assertEqual((part, acct, qos, requeue),
                             ("low", "publicgrp", "publicgrp-low-qos", True), step)


class RadiantArrayQueueTests(QueueHarness):
    def radiant(self):
        mzml = []
        for i in range(3):
            mzml.append(os.path.join(self.d, f"r{i}.mzML"))
            open(mzml[-1], "w").close()
        listing = os.path.join(self.d, "files.txt")
        with open(listing, "w") as fh:
            fh.write("\n".join(mzml) + "\n")
        fasta = os.path.join(self.d, "db.fasta")
        config = os.path.join(self.d, "default.radiantConfig")
        for p in (fasta, config):
            open(p, "w").close()
        out = os.path.join(self.d, "out")
        self.run_script([os.path.join(SCRIPTS, "radiant_parallel.py"), "--runtime", "apptainer",
                         "--image", "/x/radiant-fulcrum-2.3.3.sif", "--raw-list", listing,
                         "--fasta", fasta, "--config", config, "--out", out,
                         "--diann", "/opt/diann-2.6.1/diann-linux"])
        D = os.path.join(out, "radiant_parallel")
        return {s: queue_of(os.path.join(D, f"{s}.sbatch"))
                for s in ("step1_libpred", "step2_search", "step3_fulcrum")
                if os.path.exists(os.path.join(D, f"{s}.sbatch"))}

    def test_only_the_per_file_array_moves_to_low_when_squeue_cannot_run(self):
        self.squeue_reports(None)
        q = self.radiant()
        self.assertEqual(set(q), {"step1_libpred", "step2_search", "step3_fulcrum"}, q)
        self.assertEqual(q["step2_search"][:2], ("low", "publicgrp"))
        self.assertEqual(q["step1_libpred"][:2], ("high", "genome-center-grp"))
        self.assertEqual(q["step3_fulcrum"][:2], ("high", "genome-center-grp"))


if __name__ == "__main__":
    unittest.main()
