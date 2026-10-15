#!/usr/bin/env python3
"""
SILAC searches (estimate_params.py --labels, DIA-NN plexDIA): the parameters, and the steps they
reach.

A triple-SILAC search (2026-10-14, DIA-NN 2.7.0) had to be set up by hand with --overrides, which
recorded the study design as one person's override. Its declaration also carried
`--lib-fixed-mod SILAC`, the README's form for a label-free library: on a library predicted from a
FASTA, where `--fixed-mod SILAC,0.0,KR,label` had already labelled K and R, the precursor ids read
K(SILAC)(SILAC), and K(SILAC)(SILAC)(SILAC) once the search step applied it again (DIA-NN: "contains
more than one channel family label, this is not supported"). With no channel-normalisation flag
DIA-NN stopped and fell back to --channel-run-norm. Checked on DIA-NN 2.7.0 (HIVE srun 24808559,
2026-10-15): --fixed-mod alone labels the predicted library once.

These tests pin:
  * the masses are UniMod's, and each preset writes exactly the declaration that ran;
  * --lib-fixed-mod is never written; --channels, --original-mods and the normalisation flag are;
  * the normalisation answer is recorded, and a default is tagged as one;
  * SILAC flags passed as --overrides, a missing preset, Sage and DDA are refused;
  * a cfg that would mis-search is refused by both generators;
  * every step of the two-job route and of the 5-step chain gets the declaration, the --channels
    value reaching DIA-NN as ONE word (it holds ';'), and no step gets --lib-fixed-mod.
"""
import json
import os
import shlex
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import silac  # noqa: E402
import diann_parallel as dp  # noqa: E402
from job_env import job_env  # noqa: E402

# UniMod, read 2026-10-15 (modifications_view.php?editid1=<id>): independent of silac.py
UNIMOD = {481: 4.025107, 188: 6.020129, 259: 8.014199, 267: 10.008269}
TRIPLE = "SILAC,L,KR,0:0;SILAC,M,KR,4.025107:6.020129;SILAC,H,KR,8.014199:10.008269"


def estimate(d, *args, name="params.cfg", ok=True):
    cfg = os.path.join(d, name)
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                        "--engine", "diann", "--acquisition", "DIA", "--instrument", "timsTOF HT",
                        "--precursor-mz-range", "299.5", "1200.5", "--out", cfg, *args],
                       capture_output=True, text=True, timeout=60)
    if ok:
        assert p.returncode == 0, p.stderr
    return cfg, p


def cfg_flags(cfg):
    return {f: v for f, v in dp.cfg_groups(dp.cfg_tokens(cfg))}


class Masses(unittest.TestCase):
    def test_every_label_is_its_unimod_entry(self):
        for name, lab in silac.LABELS.items():
            if lab["unimod"] is None:
                self.assertEqual(lab["delta"], 0.0, name)
            else:
                self.assertAlmostEqual(lab["delta"], UNIMOD[lab["unimod"]], places=6, msg=name)

    def test_channel_strings_keep_every_digit(self):
        rec = silac.labels_record("silac3", "K4R6,K8R10")
        self.assertEqual(rec["diann"]["--channels"], TRIPLE)
        self.assertEqual(silac.labels_record("silac2", "K8R10")["diann"]["--channels"],
                         "SILAC,L,KR,0:0;SILAC,H,KR,8.014199:10.008269")
        self.assertEqual(silac.labels_record("silac2", "K6R6")["diann"]["--channels"],
                         "SILAC,L,KR,0:0;SILAC,H,KR,6.020129:6.020129")

    def test_presets_are_required_and_ordered(self):
        for scheme, presets in (("silac3", None), ("silac3", "K8R10"), ("silac2", "K4R6,K8R10"),
                                ("silac3", "K8R10,K4R6"), ("silac3", "K4R6,K4R6"),
                                ("silac2", "K9R9"), ("silac4", "K8R10")):
            with self.subTest(scheme=scheme, presets=presets):
                with self.assertRaises(silac.SilacError):
                    silac.labels_record(scheme, presets)


class EstimateParams(unittest.TestCase):
    def test_triple_silac_cfg_is_the_declaration_that_ran_without_lib_fixed_mod(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, p = estimate(d, "--labels", "silac3", "--label-preset", "K4R6,K8R10",
                              "--channel-norm", "spec")
            f = cfg_flags(cfg)
            self.assertEqual(f["--fixed-mod"], ["SILAC,0.0,KR,label"])
            self.assertEqual(f["--channels"], [TRIPLE])
            self.assertIn("--original-mods", f)
            self.assertIn("--channel-spec-norm", f)
            self.assertNotIn("--channel-run-norm", f)
            self.assertNotIn("--lib-fixed-mod", f)
            self.assertIn("--fasta-search", f)               # still library-free
            side = json.load(open(cfg + ".rationale.json"))
            self.assertEqual(side["overrides"], {}, "the labels are the design, not an override")
            lab = side["labels"]
            self.assertEqual([c["name"] for c in lab["channels"]], ["L", "M", "H"])
            self.assertEqual(lab["channel_norm"]["value"], "spec")
            self.assertTrue(lab["channel_norm"]["given"])
            self.assertIn("not used", side["rationale"]["--lib-fixed-mod"]["value"])
            self.assertEqual(silac.cfg_problems(dp.cfg_groups(dp.cfg_tokens(cfg))), [])

    def test_the_normalisation_default_is_tagged(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "--labels", "silac2", "--label-preset", "K8R10")
            lab = json.load(open(cfg + ".rationale.json"))["labels"]
            self.assertFalse(lab["channel_norm"]["given"])
            self.assertIn("DEFAULT, not user-confirmed", lab["channel_norm"]["source"])
            self.assertIn("DEFAULT, not user-confirmed",
                          json.load(open(cfg + ".rationale.json"))
                          ["rationale"]["--channel-spec-norm"]["source"])

    def test_turnover_writes_run_norm(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "--labels", "silac2", "--label-preset", "K8R10",
                              "--channel-norm", "run")
            f = cfg_flags(cfg)
            self.assertIn("--channel-run-norm", f)
            self.assertNotIn("--channel-spec-norm", f)
            lab = json.load(open(cfg + ".rationale.json"))["labels"]
            self.assertFalse(lab["channel_norm"]["ratio_analysis"])

    def test_refusals_leave_no_cfg(self):
        cases = [
            (["--overrides", json.dumps({"--channels": TRIPLE})], "--channels in --overrides"),
            (["--overrides", json.dumps({"--lib-fixed-mod": "SILAC"})], "--lib-fixed-mod"),
            (["--overrides", json.dumps({"--fixed-mod": "SILAC,0.0,KR,label"})], "--labels"),
            (["--labels", "silac3"], "--label-preset"),
            (["--label-preset", "K8R10"], "describe --labels"),
            (["--channel-norm", "spec"], "describe --labels"),
        ]
        for args, said in cases:
            with self.subTest(args=args), tempfile.TemporaryDirectory() as d:
                cfg = os.path.join(d, "params.cfg")
                open(cfg, "w").write("--stale\n")
                open(cfg + ".rationale.json", "w").write("{}")
                _, p = estimate(d, *args, ok=False)
                self.assertNotEqual(p.returncode, 0)
                self.assertIn(said, p.stderr)
                self.assertFalse(os.path.exists(cfg), "a stale cfg must not survive a refusal")

    def test_sage_and_dda_are_refused(self):
        with tempfile.TemporaryDirectory() as d:
            for args in (["--engine", "sage"], ["--acquisition", "DDA"]):
                cfg = os.path.join(d, "x.cfg")
                p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                                    "--engine", "diann", "--acquisition", "DIA",
                                    "--instrument", "timsTOF HT",
                                    "--precursor-mz-range", "299.5", "1200.5", *args,
                                    "--labels", "silac2", "--label-preset", "K8R10", "--out", cfg],
                                   capture_output=True, text=True, timeout=60)
                self.assertNotEqual(p.returncode, 0, args)
                self.assertIn("DIA-NN plexDIA", p.stderr)


class CfgProblems(unittest.TestCase):
    def groups(self, text):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.cfg")
            open(p, "w").write(text)
            return dp.cfg_groups(dp.cfg_tokens(p))

    BASE = ("--fasta-search --gen-spec-lib --predictor --fixed-mod SILAC,0.0,KR,label "
            f"--channels {TRIPLE} --original-mods ")

    def test_the_prot_0807_type_declaration_is_refused(self):
        # as it ran: --lib-fixed-mod in the cfg, and no normalisation flag
        probs = silac.cfg_problems(self.groups(self.BASE + "--lib-fixed-mod SILAC\n"))
        self.assertTrue(any("--lib-fixed-mod" in p for p in probs), probs)
        self.assertTrue(any("normalisation" in p for p in probs), probs)

    def test_other_problems(self):
        self.assertEqual(silac.cfg_problems(self.groups(self.BASE + "--channel-spec-norm\n")), [])
        self.assertTrue(silac.cfg_problems(self.groups(
            self.BASE + "--channel-spec-norm --channel-run-norm\n")))
        self.assertTrue(silac.cfg_problems(self.groups(
            self.BASE.replace("--original-mods ", "") + "--channel-spec-norm\n")))
        self.assertTrue(silac.cfg_problems(self.groups(
            self.BASE.replace("--fixed-mod SILAC,0.0,KR,label ", "") + "--channel-spec-norm\n")))
        self.assertTrue(silac.cfg_problems(self.groups(self.BASE + "--channel-spec-norm\n"),
                                           seed_lib=True))
        self.assertTrue(silac.cfg_problems(self.groups(
            "--fixed-mod SILAC,0.0,KR,label --original-mods --channel-spec-norm "
            "--channels SILAC,L,KR,0:0;SILAC,H,KR,8.014199:6.0;SILAC,X,KR,4:4\n")))
        # label-free: nothing to say
        self.assertEqual(silac.cfg_problems(self.groups("--fasta-search --unimod4\n")), [])

    def test_labels_read_back_from_any_cfg(self):
        rec = silac.labels_from_cfg(self.groups(self.BASE + "--lib-fixed-mod SILAC\n"))
        self.assertEqual([c["preset"] for c in rec["channels"]], ["K0R0", "K4R6", "K8R10"])
        self.assertEqual(rec["channels"][1]["K"]["unimod"], 481)
        self.assertIsNone(rec["channel_norm"])
        self.assertTrue(rec["lib_fixed_mod"])


class Routes(unittest.TestCase):
    """The flags reach the steps: the two-job route and the 5-step chain, as generated."""

    def _cfg(self, d):
        cfg, _ = estimate(d, "--labels", "silac3", "--label-preset", "K4R6,K8R10",
                          "--channel-norm", "spec")
        return cfg

    def assert_declared(self, words, where):
        self.assertIn("--channels", words, where)
        self.assertEqual(words[words.index("--channels") + 1], TRIPLE, where)
        self.assertEqual(words[words.index("--fixed-mod") + 1], "SILAC,0.0,KR,label", where)
        self.assertIn("--original-mods", words, where)
        self.assertIn("--channel-spec-norm", words, where)
        self.assertNotIn("--lib-fixed-mod", words, where)

    def test_two_job_route(self):
        import test_single_shot_mass_accuracy as ss
        with tempfile.TemporaryDirectory() as d:
            raws, fasta, _, tools, bundle = ss.SingleShotMassAccTests._setup(self, d)
            cfg = self._cfg(d)
            out = os.path.join(d, "out")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                                "--tools", tools, "--bundle", bundle, "--params", cfg,
                                "--fasta", fasta, "--out", out, "--files", *raws,
                                "--engine", "diann", "--threads", "8", "--sbatch", "job.sh",
                                "--no-notify", "--no-fran"],
                               capture_output=True, text=True, cwd=d, timeout=240,
                               env=ss.SingleShotMassAccTests._env(self, d))
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            lib = open(os.path.join(d, "job_1_lib.sh")).read()
            srch = open(os.path.join(d, "job_2_search.sh")).read()
            # job 1 reads the whole cfg (--cfg): what DIA-NN reads is the cfg itself
            self.assertIn(f"--cfg {shlex.quote(os.path.abspath(cfg))}", lib)
            self.assert_declared(dp.cfg_tokens(cfg), "the cfg job 1 reads")
            line = next(ln for ln in srch.splitlines() if " --lib " in ln)
            self.assert_declared(shlex.split(line), "job 2 (search)")

    def test_a_bad_cfg_is_refused_before_anything_is_written(self):
        import test_single_shot_mass_accuracy as ss
        with tempfile.TemporaryDirectory() as d:
            raws, fasta, _, tools, bundle = ss.SingleShotMassAccTests._setup(self, d)
            cfg = self._cfg(d)
            with open(cfg, "a") as fh:
                fh.write("--lib-fixed-mod SILAC\n")
            out = os.path.join(d, "out")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                                "--tools", tools, "--bundle", bundle, "--params", cfg,
                                "--fasta", fasta, "--out", out, "--files", *raws,
                                "--engine", "diann", "--threads", "8", "--sbatch", "job.sh"],
                               capture_output=True, text=True, cwd=d, timeout=240,
                               env=ss.SingleShotMassAccTests._env(self, d))
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("--lib-fixed-mod", p.stdout + p.stderr)
            self.assertFalse(os.path.exists(os.path.join(d, "job_2_search.sh")))

    def test_five_step_chain(self):
        from test_cfg_reader_quoting import _fake_diann, _calls
        with tempfile.TemporaryDirectory() as d:
            raws = []
            for i in range(6):
                raws.append(os.path.join(d, f"s{i}.d"))
                os.makedirs(raws[-1])
            fasta = os.path.join(d, "db.fasta")
            open(fasta, "w").write(">sp|P1|X\nPEPTIDEKR\n")
            cfg = self._cfg(d)
            fake, out = _fake_diann(d), os.path.join(d, "out")
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", fake, "--raw", *raws, "--fasta", fasta, "--out", out,
                                "--cfg", cfg, "--no-notify", "--no-fran"],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(g.returncode, 0, g.stderr)
            steps = ("step1_libpred", "step2_firstpass", "step3_assembly", "step4_finalpass",
                     "step5_report")
            for s in steps:
                text = open(os.path.join(out, s + ".sbatch")).read()
                self.assertIn(shlex.quote(TRIPLE), text, f"{s} lacks --channels")
                self.assertNotIn("--lib-fixed-mod", text, s)
                self.assertIn("--channel-spec-norm", text, s)
            # step 1 run against a fake DIA-NN: the argv it receives
            env = job_env(d, ARGV_OUT=os.path.join(d, "argv"),
                          FAKE_TOUCH=os.path.join(out, "step1.predicted.speclib"))
            s1 = subprocess.run(["bash", os.path.join(out, "step1_libpred.sbatch")], cwd=d,
                                env=env, capture_output=True, text=True)
            self.assertEqual(s1.returncode, 0, s1.stdout + s1.stderr)
            self.assert_declared(_calls(env["ARGV_OUT"])[0], "step 1 (library prediction)")

    def test_the_chain_refuses_a_seed_library(self):
        with tempfile.TemporaryDirectory() as d:
            raws = []
            for i in range(6):
                raws.append(os.path.join(d, f"s{i}.d"))
                os.makedirs(raws[-1])
            fasta = os.path.join(d, "db.fasta")
            open(fasta, "w").write(">sp|P1|X\nPEPTIDEKR\n")
            seed = os.path.join(d, "seed.speclib")
            open(seed, "w").write("x")
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", "/bin/true", "--raw", *raws, "--fasta", fasta,
                                "--out", os.path.join(d, "out"), "--cfg", self._cfg(d),
                                "--seed-lib", seed], capture_output=True, text=True, timeout=120)
            self.assertNotEqual(g.returncode, 0)
            self.assertIn("seed library", g.stderr)


class MethodsSentence(unittest.TestCase):
    def test_masses_unimod_and_normalisation_answer(self):
        s = silac.methods_search_sentence(silac.labels_record("silac3", "K4R6,K8R10", "spec"))
        for want in ("plexDIA", "Derks", "Lys +4.025107 Da (2H(4), UniMod:481)",
                     "Arg +10.008269 Da (13C(6)15N(4), UniMod:267)", "--channel-spec-norm",
                     "L (unlabelled Lys and Arg)"):
            self.assertIn(want, s)
        self.assertNotIn("DEFAULT", s)
        d = silac.methods_search_sentence(silac.labels_record("silac2", "K8R10"))
        self.assertIn("DEFAULT, not user-confirmed", d)
        self.assertEqual(silac.methods_search_sentence(None), "")


if __name__ == "__main__":
    unittest.main()
