#!/usr/bin/env python3
"""
A DIA-NN cfg is read ONCE, by one rule, and its values reach DIA-NN unchanged -- on both routes.

The cfg's flags are spliced into bash: every step of the 5-step chain, and the single-shot
library-free search. So two things have to hold for any cfg the parallel gate approves:

  1. the reader splits words the way bash does (diann_parallel.cfg_tokens), so the gate, the
     step flags, params.base.cfg and run_search all see the same flags; and
  2. the emitter (diann_parallel._shield / bash_flags) turns each word back into ONE bash word
     whose value is exactly that word.

What broke, each confirmed:
  * `--var-mod "Phospho(STY),79.966331,STY"` lost its quotes in the reader and was re-emitted
    bare: every step died on "syntax error near unexpected token `('" -- after the gate had
    approved the cfg. Phospho is a common modification.
  * shlex ended a word at ANY `#` (`/data/run#1` -> `/data/run`); bash keeps it.
  * `x{1,2}` and `~/libs` were emitted bare, so bash rewrote them.
  * run_search's single-shot search used `cfg_txt.split()`: a `# comment` commented out
    --f/--fasta/--out, `--cut K*,R*` went in bare, and ensure_xic counted a commented-out
    `# --xic 10` as present -- while the chain's reader did not, and step 4 extracted nothing.

The cfgs estimate_params.py writes must come out exactly as before (they are what actually runs
day to day): the chain flags byte-identical, the single-shot command byte-identical except the
one quoted glob, and the argv DIA-NN receives identical.
"""
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
from job_env import job_env  # noqa: E402  (env for running job scripts)

import diann_parallel as dp  # noqa: E402
import run_search  # noqa: E402


def _write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def _read(path):
    with open(path) as fh:
        return fh.read()


def _bash_words(fragment, cwd=None):
    """The words bash makes of `fragment` as command arguments, with globbing and brace
    expansion OFF -- i.e. pure word splitting and quote removal. A START sentinel keeps an empty
    fragment distinguishable from one empty word."""
    p = subprocess.run(["bash", "-c", f"set -f +B\nprintf '%s\\0' START {fragment}"],
                       capture_output=True, text=True, cwd=cwd)
    assert p.returncode == 0, p.stderr
    words = p.stdout.split("\0")[:-1]
    assert words[0] == "START"
    return words[1:]


def _bash_argv(fragment, cwd=None, env=None):
    """The words bash hands a command for `fragment`, with ALL expansions on (as in a step)."""
    p = subprocess.run(["bash", "-c", f"printf '%s\\0' START {fragment}"],
                       capture_output=True, text=True, cwd=cwd, env=env)
    assert p.returncode == 0, p.stderr
    return p.stdout.split("\0")[1:-1]


class ReaderMatchesBashTests(unittest.TestCase):

    LINES = [
        "--qvalue 0.01   # 1% FDR",
        "--path /data/run#1 --x 2",
        '--var-mod "Phospho(STY),79.966331,STY"',
        "--a 'single quoted # not a comment' --b \"double \\\"esc\\\" \\\\ back\"",
        "--c back\\ slash\\ space",
        '--d ""',
        "--e a#b #real comment --not-a-flag",
        "#whole line comment",
        "--f 'it'\"'\"'s'",
        "--g \"keep \\n literal\"",
        "\t--h\t7  ",
    ]

    def test_the_reader_splits_words_exactly_as_bash_does(self):
        with tempfile.TemporaryDirectory() as d:
            for line in self.LINES:
                with self.subTest(line=line):
                    cfg = os.path.join(d, "p.cfg")
                    _write(cfg, line + "\n")
                    self.assertEqual(dp.cfg_tokens(cfg), _bash_words(line))

    def test_crlf_line_endings_do_not_leak_into_values(self):
        """A deliberate difference from bash: a cfg saved on Windows must not turn
        `--window 7` into `7\\r` (and fail as "not a positive integer")."""
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w", newline="") as fh:
                fh.write("--mass-acc 15\r\n--mass-acc-ms1 15\r\n--window 7\r\n")
            self.assertEqual(dp.cfg_tokens(cfg)[-2:], ["--window", "7"])
            self.assertEqual(dp.parallel_safe(cfg)["code"], "pinned")

    def test_an_unclosed_quote_is_an_error_not_a_guess(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            _write(cfg, '--var-mod "UniMod:35,15.994915,M\n')
            with self.assertRaises(dp.CfgError) as cm:
                dp.cfg_tokens(cfg)
            self.assertEqual(cm.exception.code, "cfg_unparseable")


class EmitterRoundTripTests(unittest.TestCase):

    NASTY = ["Phospho(STY),79.966331,STY", "x{1,2}", "{a,b}", "~", "~/libs", "a:~/x",
             "/data/run#1", "#lead", "a;b", "a&b", "a|b", "a<b>c", "it's", 'say "hi"',
             "back\\slash", "`id`", "$(echo pwned)", "$", "$1", "!bang", "a b", "tab\tx",
             "*", "?", "[ab]", "K*,R*", "UniMod:1,42.010565,*n", "", "=eq", "%pct", "^caret",
             "100%", "-3", "a\nb"]

    def test_every_token_reaches_bash_as_one_unchanged_word(self):
        """With EVERY expansion on and a file in the cwd that each glob would match."""
        with tempfile.TemporaryDirectory() as d:
            for name in ("Kfoo,Rbar", "UniMod:1,42.010565,Xn", "a", "b", "x"):
                _write(os.path.join(d, name), "")
            for tok in self.NASTY:
                with self.subTest(token=tok):
                    self.assertEqual(_bash_argv(dp._shield(tok), cwd=d), [tok])

    def test_estimate_params_values_stay_bare(self):
        for tok in ("--qvalue", "0.01", "UniMod:35,15.994915,M", "Cont_", "--mass-acc-ms1",
                    "299", "4.0"):
            self.assertEqual(dp._shield(tok), tok)
        self.assertEqual(dp._shield("K*,R*"), "'K*,R*'")        # as since 2.4.2

    def test_only_name_style_variables_expand(self):
        env = dict(os.environ, HOME="/home/test-user")
        self.assertEqual(_bash_argv(dp._shield("$HOME/libs"), env=env), ["/home/test-user/libs"])
        self.assertEqual(_bash_argv(dp._shield("${HOME}/x (y)"), env=env), ["/home/test-user/x (y)"])
        self.assertEqual(_bash_argv(dp._shield("$HOME/`id`$(id)"), env=env),
                         ["/home/test-user/`id`$(id)"])


FAKE_DIANN_ARGV = """#!/bin/bash
printf '%s\\0' "$@" >> "$ARGV_OUT"
printf '\\n\\0' >> "$ARGV_OUT"
[ -n "$FAKE_TOUCH" ] && echo lib > "$FAKE_TOUCH"
echo "Scan window radius set to 7"
"""


def _fake_diann(d):
    fake = os.path.join(d, "fake-diann")
    _write(fake, FAKE_DIANN_ARGV)
    os.chmod(fake, 0o755)
    return fake


def _calls(argv_out):
    """Each fake-DIA-NN invocation's argv, split on the newline record marker."""
    calls, cur = [], []
    for w in _read(argv_out).split("\0")[:-1]:
        if w == "\n":
            calls.append(cur)
            cur = []
        else:
            cur.append(w)
    return calls


class PhosphoCfgRunsTests(unittest.TestCase):
    """Round 2, item 1, end to end: generate the chain for a cfg with every awkward value, check
    every script parses, then RUN step 1 and step 1b against a fake DIA-NN that records argv."""

    CFG = ("--qvalue 0.01   # 1% FDR\n"
           "--mass-acc 15 --mass-acc-ms1 15\n"
           '--var-mod "Phospho(STY),79.966331,STY"\n'
           '--note "x{1,2}" --home ~/libs --path /data/run#1\n'
           "--cut K*,R*\n")
    WANT = ["--qvalue", "0.01", "--mass-acc", "15", "--mass-acc-ms1", "15",
            "--var-mod", "Phospho(STY),79.966331,STY", "--note", "x{1,2}", "--home", "~/libs",
            "--path", "/data/run#1", "--cut", "K*,R*"]

    def test_steps_parse_and_dianns_argv_is_the_cfg(self):
        with tempfile.TemporaryDirectory() as d:
            raws = []
            for i in range(6):
                raws.append(os.path.join(d, f"s{i}.d"))
                os.makedirs(raws[-1])
            fasta = os.path.join(d, "db.fasta")
            _write(fasta, ">sp|P1|X\nPEPTIDER\n")
            cfg = os.path.join(d, "p.cfg")
            _write(cfg, self.CFG)
            for name in ("Kfoo,Rbar", "x1", "x2"):           # things bash would expand into
                _write(os.path.join(d, name), "")
            fake, out = _fake_diann(d), os.path.join(d, "out")
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", fake, "--raw", *raws, "--fasta", fasta, "--out", out,
                                "--cfg", cfg], capture_output=True, text=True)
            self.assertEqual(g.returncode, 0, g.stderr)
            self.assertTrue(dp.parallel_safe(cfg)["ok"])
            for f in sorted(os.listdir(out)):
                if f.endswith((".sbatch", ".sh")):
                    n = subprocess.run(["bash", "-n", os.path.join(out, f)],
                                       capture_output=True, text=True)
                    self.assertEqual(n.returncode, 0, f"{f}: {n.stderr}")

            env = job_env(d, ARGV_OUT=os.path.join(d, "argv"),
                          FAKE_TOUCH=os.path.join(out, "step1.predicted.speclib"))
            s1 = subprocess.run(["bash", os.path.join(out, "step1_libpred.sbatch")], cwd=d,
                                env=env, capture_output=True, text=True)
            self.assertEqual(s1.returncode, 0, s1.stdout + s1.stderr)
            env.pop("FAKE_TOUCH")
            s1b = subprocess.run(["bash", os.path.join(out, "step1b_window.sbatch")], cwd=d,
                                 env=env, capture_output=True, text=True, timeout=120)
            self.assertEqual(s1b.returncode, 0, s1b.stdout + s1b.stderr)
            step1, *probes = _calls(env["ARGV_OUT"])
            self.assertEqual(len(probes), dp.PROBE_CANDIDATES)     # representative runs
            for label, argv in [("step 1", step1)] + [("step 1b probe", pr) for pr in probes]:
                k = argv.index("--threads") + 2
                self.assertEqual(argv[k:], self.WANT, f"{label} handed DIA-NN different flags")
            # and the resolved cfg reads back to the same values plus the measured window
            self.assertEqual(dp.cfg_tokens(os.path.join(out, "params.resolved.cfg")),
                             self.WANT + ["--window", "7"])


def _estimate(d, name, *args):
    cfg = os.path.join(d, f"{name}.cfg")
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                        "--engine", "diann", "--out", cfg, *args], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return cfg


def _estimate_params_cfgs(d):
    meta = os.path.join(d, "db.fasta.meta.json")
    _write(meta, json.dumps({"diann_cont_quant_exclude": "Cont_"}))
    return {
        "timsTOF": _estimate(d, "tims", "--acquisition", "DIA", "--instrument", "timsTOF HT",
                             "--precursor-mz-range", "299.5", "1200.5", "--fasta-meta", meta),
        "Astral": _estimate(d, "astral", "--acquisition", "DIA", "--instrument", "Orbitrap Astral"),
        "Orbitrap": _estimate(d, "orbi", "--acquisition", "DIA", "--instrument", "Exploris 480",
                              "--ms1-resolution", "120000", "--ms2-resolution", "30000",
                              "--var-mods", "ox"),
    }


class SameAsOriginForEstimateParamsTests(unittest.TestCase):
    """What worked must not change. Origin's rules are reproduced verbatim from df92b95."""

    @staticmethod
    def origin_chain_flags(cfg, drop=()):
        globby = re.compile(r"[*?\[\]]")
        out, strip = [], tuple(dp.STRIP) + tuple(drop)
        for raw in open(cfg):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if any(line == s or line.startswith(s + " ") for s in strip):
                continue
            out.append(" ".join(shlex.quote(t) if globby.search(t) else t for t in line.split()))
        return " ".join(out)

    ORIGIN_SINGLE_SHOT_STRIP = ("--fasta-search", "--gen-spec-lib", "--predictor", "--reanalyse",
                                "--matrices", "--rt-profiling")

    @staticmethod
    def origin_single_shot(cmd, params, files, fasta, out, threads, dda="",
                           strip=ORIGIN_SINGLE_SHOT_STRIP):
        report = os.path.join(out, "report.parquet")
        f_args = " ".join(f"--f {shlex.quote(f)}" for f in files)
        lib = os.path.join(out, "diann_lib")
        search_cfg = " ".join(t for t in open(params).read().split() if t not in strip)
        return (f"{cmd} {search_cfg} {f_args} --fasta {shlex.quote(fasta)} "
                f"--lib {shlex.quote(lib)}.predicted.speclib --reanalyse --matrices "
                f"--out {shlex.quote(report)} --threads {threads}{dda}")

    def test_chain_flags_are_byte_identical(self):
        with tempfile.TemporaryDirectory() as d:
            for name, cfg in _estimate_params_cfgs(d).items():
                with self.subTest(name):
                    self.assertEqual(dp.read_cfg_flags(cfg), self.origin_chain_flags(cfg))
                    self.assertEqual(dp.read_cfg_flags(cfg, drop=("--window",)),
                                     self.origin_chain_flags(cfg, drop=("--window",)))

    def _single_shot(self, d, cfg, fake):
        files = [os.path.join(d, "a b.d"), os.path.join(d, "s1.d")]
        out, fasta = os.path.join(d, "out"), os.path.join(d, "db.fasta")
        job = os.path.join(d, "job.sh")
        run_search.run_diann(fake, cfg, files, fasta, out, 8, job, acquisition="DIA")
        # The DIA-NN line, not the job's last line: the job also clears and checks its report
        # around it (tests/test_single_shot_sbatch.py).
        lines = [l for l in _read(os.path.join(d, "job_2_search.sh")).splitlines()
                 if l.startswith(fake + " ")]
        self.assertEqual(len(lines), 1, lines)
        old = self.origin_single_shot(fake, cfg, files, fasta, out, 8)
        # The search job now KEEPS --rt-profiling: with --reanalyse it builds the empirical
        # library that flag configures (run_search.SINGLE_SHOT_SEARCH_STRIP says why).
        kept = self.origin_single_shot(fake, cfg, files, fasta, out, 8, strip=tuple(
            f for f in self.ORIGIN_SINGLE_SHOT_STRIP if f != "--rt-profiling"))
        return lines[0], old, kept

    def test_single_shot_command_is_byte_identical_but_for_the_quoted_glob(self):
        """(a) the intended differences are `K*,R*` -> `'K*,R*'`, exactly as the chain has
        quoted it since 2.4.2, and `--rt-profiling` kept in place. (b) bash hands DIA-NN the
        same argv plus that one flag -- and (c) with a file matching the glob in the working
        directory, origin's bare form is rewritten while the new one is not."""
        with tempfile.TemporaryDirectory() as d:
            fake = _fake_diann(d)
            for name, cfg in _estimate_params_cfgs(d).items():
                with self.subTest(name):
                    new, old, kept = self._single_shot(d, cfg, fake)
                    # ...and one more intended difference since 2.5.1: the search names its own
                    # --temp (<out>/quant) when the cfg has none, so DIA-NN does not write .quant
                    # files beside the raw data (run_diann, "--temp, always").
                    new, n_temp = re.subn(r" --temp \S+/quant$", "", new)
                    self.assertEqual(n_temp, 1, "the single-shot search names no --temp")
                    self.assertEqual(old.count(" K*,R* "), 1)
                    self.assertNotEqual(kept, old, "estimate_params wrote no --rt-profiling")
                    self.assertEqual(new, kept.replace(" K*,R* ", " 'K*,R*' "))

                    clean = os.path.join(d, f"clean-{name}")
                    os.makedirs(clean)
                    got = _bash_argv(new.split(" ", 1)[1], cwd=clean)
                    self.assertEqual(got.count("--rt-profiling"), 1)
                    self.assertEqual([t for t in got if t != "--rt-profiling"],
                                     _bash_argv(old.split(" ", 1)[1], cwd=clean))

                    trap = os.path.join(d, f"trap-{name}")
                    os.makedirs(trap)
                    _write(os.path.join(trap, "Kfoo,Rbar"), "")
                    self.assertIn("Kfoo,Rbar", _bash_argv(old.split(" ", 1)[1], cwd=trap),
                                  "origin's bare glob should have been rewritten here")
                    got = _bash_argv(new.split(" ", 1)[1], cwd=trap)
                    self.assertIn("K*,R*", got)
                    self.assertNotIn("Kfoo,Rbar", got)

    def test_a_comment_no_longer_swallows_the_single_shot_search(self):
        with tempfile.TemporaryDirectory() as d:
            fake = _fake_diann(d)
            cfg = os.path.join(d, "c.cfg")
            _write(cfg, "--fasta-search --gen-spec-lib   # library-free\n--qvalue 0.01 # 1%\n"
                        "--mass-acc 15\n")
            new, _, _ = self._single_shot(d, cfg, fake)
            argv = _bash_argv(new.split(" ", 1)[1], cwd=d)
            for flag in ("--qvalue", "--mass-acc", "--f", "--fasta", "--lib", "--out", "--threads"):
                self.assertIn(flag, argv, f"{flag} was commented out of the single-shot search")
            self.assertNotIn("#", " ".join(argv))


class EnsureXicReadsLikeTheChainTests(unittest.TestCase):

    def test_a_commented_out_xic_is_not_counted_as_present(self):
        """split() saw `--xic` inside `# --xic 10`, added nothing, and the chain -- reading the
        same cfg properly -- found no --xic: step 4 extracted no chromatograms, silently."""
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            _write(cfg, "--qvalue 0.01\n# --xic 10\n")
            used = run_search.ensure_xic(cfg, os.path.join(d, "out"))
            self.assertNotEqual(used, cfg)
            self.assertIn("--xic 10", dp.xic_flag(used))

    def test_estimate_params_cfgs_are_still_left_alone(self):
        with tempfile.TemporaryDirectory() as d:
            for name, cfg in _estimate_params_cfgs(d).items():
                with self.subTest(name):
                    self.assertEqual(run_search.ensure_xic(cfg, os.path.join(d, "out")), cfg)


if __name__ == "__main__":
    unittest.main(verbosity=2)
