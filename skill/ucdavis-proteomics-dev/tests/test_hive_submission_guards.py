#!/usr/bin/env python3
"""
The two things that cost a HIVE submission cycle.

Both failures are cheap to prevent and expensive to hit: you queue, the job starts, dies in
seconds, and you read a log and resubmit. Real users lost cycles to each.

  1. `--temp` must pre-exist. DIA-NN exits with "cannot find the temp folder" and will NOT
     create it — before doing any work.
  2. The queue must be the one this account can actually submit to. A hardcoded facility queue
     is REJECTED for an account outside the group; a hardcoded preemptible queue silently
     exposes a long job to preemption.
"""
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import run_search  # noqa: E402


class TempDirTests(unittest.TestCase):
    def test_a_temp_dir_named_in_the_cfg_is_created(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "search_out")
            os.makedirs(out)
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write(f"--qvalue 0.01\n--temp {os.path.join(d, 'scratch')}\n")
            made = run_search.ensure_temp_dirs(cfg, out)
            self.assertTrue(os.path.isdir(os.path.join(d, "scratch")))
            self.assertEqual(made, [os.path.join(d, "scratch")])

    def test_a_relative_temp_resolves_against_the_output_dir(self):
        """DIA-NN runs with the output directory as its working context, so a bare `--temp tmp`
        means <out>/tmp. Creating it beside the cfg instead would leave DIA-NN's own path missing
        and the guard would look like it had worked."""
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "search_out")
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("--temp diann_tmp\n")
            run_search.ensure_temp_dirs(cfg, out)
            self.assertTrue(os.path.isdir(os.path.join(out, "diann_tmp")))

    def test_no_temp_in_the_cfg_creates_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n")
            self.assertEqual(run_search.ensure_temp_dirs(cfg, os.path.join(d, "o")), [])

    def test_it_is_idempotent_and_survives_a_missing_cfg(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write(f"--temp {os.path.join(d, 'x')}\n")
            run_search.ensure_temp_dirs(cfg, d)
            run_search.ensure_temp_dirs(cfg, d)          # second call must not raise
            self.assertEqual(run_search.ensure_temp_dirs("/no/such.cfg", d), [])


class ParallelChainTempTests(unittest.TestCase):
    """Every step of the 5-step chain makes its own temp dir.

    submit.sh already did it, but references/watcher.md tells the orchestrator to resubmit
    INDIVIDUAL steps after a failure ("resubmit steps 3→4→5 fresh"), and `sbatch
    step4_finalpass.sbatch` never runs submit.sh.
    """

    def _chain(self, d):
        raws = []
        for i in range(3):
            r = os.path.join(d, f"f{i}.d")
            os.makedirs(r)
            raws.append(r)
        fasta = os.path.join(d, "db.fasta")
        with open(fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDER\n")
        cfg = os.path.join(d, "p.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n--window 7\n")
        out = os.path.join(d, "out")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                            "--diann", "/bin/true", "--raw", *raws, "--fasta", fasta,
                            "--out", out, "--cfg", cfg],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        return out

    def test_every_step_that_uses_temp_also_creates_it(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._chain(d)
            for step in ("step2_firstpass", "step3_assembly", "step4_finalpass", "step5_report"):
                path = os.path.join(out, f"{step}.sbatch")
                self.assertTrue(os.path.exists(path), step)
                body = open(path).read()
                # the guard line itself mentions --temp in its comment; look only at the
                # DIA-NN invocation lines
                temp = re.findall(r"--temp\s+(\S+)", body.replace("\n", " "))
                temp = [t for t in temp if t.startswith("/")]
                self.assertTrue(temp, f"{step} has no --temp path")
                # the dir each step hands DIA-NN must be mkdir'd by that same step
                for which in set(temp):
                    self.assertIn(f"mkdir -p {which}", body,
                                  f"{step} passes --temp {which} but never creates it")

    def test_submit_sh_still_makes_them_too(self):
        """Belt and braces: the per-step guard is the fix, but submit.sh creating them up front
        means the first run never depends on step ordering."""
        with tempfile.TemporaryDirectory() as d:
            body = open(os.path.join(self._chain(d), "submit.sh")).read()
            self.assertIn("mkdir -p", body)
            self.assertIn("quant_step2", body)
            self.assertIn("quant_step4", body)


class CfgSplicingTests(unittest.TestCase):
    """cfg values reach DIA-NN through a bash command line, so bash gets a say in them.

    DIA-NN's own recommended values contain glob characters: `--cut K*,R*` IS the trypsin rule,
    and an N-terminal variable mod is spelled `--var-mod UniMod:1,42.010565,*n`. Unquoted, those
    are subject to pathname expansion. With nothing matching they pass through untouched, which
    is why this went unnoticed -- but one file in the working directory matching the pattern
    rewrites them silently, and the search then reports success with a changed digest rule.
    """

    def _flags(self, lines, files=()):
        import diann_parallel
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("\n".join(lines) + "\n")
            for f in files:
                open(os.path.join(d, f), "w").close()
            flags = diann_parallel.read_cfg_flags(cfg)
            # what bash actually hands the engine, evaluated in that directory
            argv = subprocess.run(["bash", "-c", f"cd {d} && printf '%s\\n' {flags}"],
                                  capture_output=True, text=True).stdout.split("\n")
            return flags, [x for x in argv if x]

    def test_a_glob_value_survives_a_matching_file_in_the_working_dir(self):
        _, argv = self._flags(["--qvalue 0.01", "--cut K*,R*"], files=["Kfoo,Rbar"])
        self.assertIn("K*,R*", argv, "bash rewrote the digest rule from a filename")
        self.assertNotIn("Kfoo,Rbar", argv)

    def test_an_nterm_variable_mod_survives_too(self):
        mod = "UniMod:1,42.010565,*n"
        _, argv = self._flags([f"--var-mod {mod}"], files=["UniMod:1,42.010565,Xn"])
        self.assertIn(mod, argv)

    def test_shell_variables_still_expand(self):
        """Only globby tokens are quoted. Quoting everything would disable $VAR and $(...) for
        every existing cfg -- a bigger behaviour change than the bug being fixed."""
        _, argv = self._flags(["--qvalue 0.01", "--lib-dir $HOME/libs"])
        self.assertIn(os.path.expanduser("~/libs"), argv)

    def test_ordinary_values_are_not_quoted_at_all(self):
        """The generated line stays byte-identical for a cfg with nothing at risk, so this
        cannot change what any existing search runs."""
        flags, _ = self._flags(["--qvalue 0.01", "--mass-acc 15"])
        self.assertEqual(flags, "--qvalue 0.01 --mass-acc 15")


class VersionBumpTests(unittest.TestCase):
    """The version lives in three places and CI fails the run if any disagrees. Two consecutive
    PRs bumped plugin.json alone; nothing about the third location is discoverable from the file
    you are editing."""

    def _bump(self, *args):
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, "bump_version.py"), *args],
                              capture_output=True, text=True)

    def test_check_reports_agreement_and_exits_zero(self):
        p = self._bump("--check")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("agree: True", p.stdout)

    def test_a_bad_version_is_refused(self):
        self.assertNotEqual(self._bump("2.4").returncode, 0)
        self.assertNotEqual(self._bump("banana").returncode, 0)


class QueueTests(unittest.TestCase):
    def test_an_explicit_queue_is_always_honoured(self):
        self.assertEqual(run_search.slurm_queue("high", "genome-center-grp", "x"),
                         ("high", "genome-center-grp", "x"))

    def test_nothing_is_invented_when_slurm_cannot_be_asked(self):
        """Emitting a guessed queue is worse than emitting none: on HIVE the cluster default
        partition is `high`, which is exactly what a non-facility account cannot use. With no
        answer, omit the lines and let SLURM apply its own."""
        part, acct, qos = run_search.slurm_queue()      # no sacctmgr on a laptop
        if part is None:
            self.assertEqual((part, acct, qos), (None, None, None))
        else:                                            # on a real cluster: a usable pair
            self.assertIn((part, acct), [("high", "genome-center-grp"), ("low", "publicgrp")])

    def test_every_sbatch_emitter_derives_its_queue_from_one_place(self):
        """diatracer_parallel.py used to hardcode --partition low with no --account and no
        --qos: facility members were put on the preemptible queue for a ~20-min-per-file
        conversion, and publicgrp/low needs publicgrp-low-qos to be accepted at all."""
        for name in ("run_search.py", "diann_parallel.py", "radiant_parallel.py",
                     "diatracer_parallel.py"):
            src = open(os.path.join(SCRIPTS, name)).read()
            if "#SBATCH" not in src:
                continue
            self.assertIn("slurm_queue", src, f"{name} emits #SBATCH without queue detection")

    def test_no_sbatch_emitter_demotes_to_the_preemptible_queue_in_silence(self):
        """A fallback to publicgrp/low is a real demotion for a facility member — a multi-hour
        search newly exposed to preemption. It happened silently when queue detection raised,
        which reads as "the skill just uses low". Every fallback has to say so on stderr."""
        for name in ("diann_parallel.py", "radiant_parallel.py", "diatracer_parallel.py"):
            src = open(os.path.join(SCRIPTS, name)).read()
            i = src.index("slurm_queue")
            tail = src[i:i + 2000]
            self.assertIn("except Exception", tail, name)
            block = tail[tail.index("except Exception"):]
            self.assertTrue("stderr.write" in block[:900],
                            f"{name} falls back to a different queue without saying so")

    def test_high_needs_no_explicit_qos(self):
        """Measured on HIVE 2026-08-25: a job submitted with --partition high --account
        genome-center-grp and NO --qos is accepted and SLURM assigns genome-center-grp-high-qos.
        So omitting it is correct, not an oversight — this test exists so nobody "fixes" the
        generator by hardcoding a QOS that a non-facility account cannot use."""
        src = open(os.path.join(SCRIPTS, "diann_parallel.py")).read()
        self.assertNotIn('qos = "genome-center-grp-high-qos"', src)
        self.assertIn("publicgrp-low-qos", src)      # low DOES need its qos named

    def test_diatracer_emits_qos_and_requeues_only_on_the_preemptible_queue(self):
        src = open(os.path.join(SCRIPTS, "diatracer_parallel.py")).read()
        self.assertIn("--qos=", src)
        # --requeue is meaningful only where preemption happens; conditioning it on `low`
        # keeps a non-preemptible submission honest about what it is.
        self.assertIn('if a.partition == "low"', src)
        self.assertNotIn('ap.add_argument("--partition", default="low")', src)


if __name__ == "__main__":
    unittest.main(verbosity=2)
