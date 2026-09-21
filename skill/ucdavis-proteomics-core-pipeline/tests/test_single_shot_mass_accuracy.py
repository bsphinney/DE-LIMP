#!/usr/bin/env python3
"""
A single-shot DIA-NN search measures an undocumented Orbitrap mass accuracy too.

Measuring only inside the 5-step chain left every machine without SLURM -- a laptop, a
workstation, a non-SLURM cluster -- with DIA-NN's first-run auto mode for the same data:
run_search.parallel_decision() returns "no SLURM here" before it looks at the cfg, so those
cohorts go single-shot at ANY size. DIA-NN 2.7.0 labels that mode "use this mode for preliminary
analyses only", the result depends on which file sorts first, and search_provenance.json never
records the value DIA-NN chose. SKILL.md golden rule 7: HIVE is a fast path, never a
requirement. probe_window.py has no SLURM dependency, so the single-shot search runs it too:
after the library job and before the search, on the same representative runs, pinning what it
measures (and the documented level as documented) for the search.

run_search.py is run as the orchestrator runs it, against a fake DIA-NN that predicts a library,
replays the captured HIVE logs for a one-run probe, and writes a report for the search.
"""
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
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import diann_parallel as dp  # noqa: E402
import run_search  # noqa: E402
import test_step1b_mass_accuracy_probe as fx  # noqa: E402  (captured logs and the probe fake)

# library job (--out-lib) -> write the predicted library; one --f -> the probe fake, replaying the
# run's captured log; several --f -> the search: write --out.
DISPATCH = r"""#!/bin/bash
nf=0; outlib=""; out=""; prev=""
for a in "$@"; do
  [ "$a" = --f ] && nf=$((nf + 1))
  [ "$prev" = --out-lib ] && outlib="$a"
  [ "$prev" = --out ] && out="$a"
  prev="$a"
done
if [ -n "$outlib" ]; then
  [ -n "${FAKE_ARGV_LOG:-}" ] && echo "LIB $*" >> "$FAKE_ARGV_LOG"
  echo lib > "$outlib.predicted.speclib"; exit 0
fi
if [ "$nf" -gt 1 ]; then
  [ -n "${FAKE_ARGV_LOG:-}" ] && echo "SEARCH $*" >> "$FAKE_ARGV_LOG"
  echo report > "$out"; exit 0
fi
FAKE_ARGV_LOG="${FAKE_ARGV_LOG:+$FAKE_ARGV_LOG.probe}" exec "$(dirname "$0")/probe_diann" "$@"
"""


class SingleShotMassAccTests(unittest.TestCase):

    def _setup(self, d, ms2_res=15000, runs=None):
        raws = runs or fx._cohort(d)
        fasta, _ = fx._fasta_lib(d)
        cfg = os.path.join(d, "params.cfg")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "diann", "--acquisition", "DIA",
                            "--instrument", "Orbitrap Exploris 480",
                            "--ms1-resolution", "120000", "--ms2-resolution", str(ms2_res),
                            "--precursor-mz-range", "357", "1105", "--out", cfg],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        fx._exe(os.path.join(d, "probe_diann"), fx.FAKE_DIANN)
        diann = fx._exe(os.path.join(d, "diann"), DISPATCH)
        tools = os.path.join(d, "tools.json")
        with open(tools, "w") as fh:
            json.dump({"diann": diann}, fh)
        bundle = os.path.join(d, "bundle.json")
        with open(bundle, "w") as fh:
            json.dump({"acquisition": "DIA", "engine": {"name": "diann"}}, fh)
        return raws, fasta, cfg, tools, bundle

    def _env(self, d):
        env = {k: v for k, v in os.environ.items() if k != "DOTNET_ROOT"}
        root = os.path.join(d, "dotnet8")
        env.update(PROTEOMICS_DOTNET_DIR=root if os.path.isdir(root) else fx._fake_dotnet_root(d),
                   FAKE_ARGV_LOG=os.path.join(d, "argv.txt"), FAKE_SEARCH_SLEEP="0",
                   # no sbatch: a machine without SLURM (python3 stays reachable)
                   PATH=os.pathsep.join([os.path.dirname(sys.executable), "/usr/bin", "/bin"]))
        return env

    def _run_search(self, d, raws, fasta, cfg, tools, bundle, *more):
        out = os.path.join(d, "out")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                            "--tools", tools, "--bundle", bundle, "--params", cfg,
                            "--fasta", fasta, "--out", out, "--files", *raws,
                            "--engine", "diann", "--threads", "8", *more],
                           capture_output=True, text=True, env=self._env(d), timeout=240)
        return p, out

    def _calls(self, d, kind):
        path = os.path.join(d, "argv.txt")
        return [ln for ln in (open(path).read().splitlines() if os.path.exists(path) else [])
                if ln.startswith(kind + " ")]

    def test_sbatch_single_shot_measures_before_the_search_and_pins_it(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = self._run_search(d, *self._setup(d), "--sbatch", os.path.join(d, "job.sh"))
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            lib_sh = open(os.path.join(d, "job_1_lib.sh")).read()
            srch_sh = open(os.path.join(d, "job_2_search.sh")).read()
            massacc = os.path.join(out, "massacc.txt")
            resolved = os.path.join(out, "params.resolved.cfg")
            self.assertNotIn("probe_window.py", lib_sh)
            self.assertIn("probe_window.py", srch_sh)
            self.assertIn("--measure mass-acc", srch_sh)
            self.assertIn("--ms1-ppm 7", srch_sh)
            # the cfg flags reach the probe's DIA-NN as bash words after `--`, as in step 1b
            self.assertRegex(srch_sh, r"probe_window\.py.* -- --qvalue 0\.01 ")
            self.assertIn(f"$(cat {massacc})", srch_sh)
            self.assertLess(srch_sh.index("probe_window.py"), srch_sh.index(f"$(cat {massacc})"))
            self.assertFalse(os.path.exists(resolved), "a resolved cfg exists before anything "
                                                       "was measured")
            prov = json.load(open(os.path.join(out, "search_provenance.json")))
            self.assertEqual(prov["resolved_params_file"], resolved)
            self.assertEqual(prov["resolved_params_produced"], "runtime")

            env = self._env(d)
            for script in ("job_1_lib.sh", "job_2_search.sh"):
                r = subprocess.run(["bash", os.path.join(d, script)], capture_output=True,
                                   text=True, env=env, timeout=240)
                self.assertEqual(r.returncode, 0, script + "\n" + r.stdout + r.stderr)
            self.assertEqual(open(massacc).read().strip(), "--mass-acc 20 --mass-acc-ms1 7")
            # the measurement is not erased by the floor: both numbers survive in the evidence
            ev = json.load(open(os.path.join(out, "mass_acc.json")))["mass_acc"]
            self.assertEqual(ev["measured_ms2_ppm"], 14.0)
            self.assertEqual(ev["pinned_ms2_ppm"], 20.0)
            self.assertTrue(ev["floored"]["ms2_ppm"])
            self.assertEqual(ev["sop_floor"], {"ms1_ppm": 7.0, "ms2_ppm": 20.0})
            probes = open(os.path.join(d, "argv.txt.probe")).read().splitlines()
            self.assertEqual(len(probes), 3, "the representative runs, one DIA-NN each")
            for call in probes:
                self.assertNotIn("--mass-acc", call, "a probe ran with mass accuracy pinned")
            search = self._calls(d, "SEARCH")
            self.assertEqual(len(search), 1)
            self.assertIn("--mass-acc 20 --mass-acc-ms1 7", search[0])
            resolved = open(os.path.join(out, "params.resolved.cfg")).read()
            self.assertEqual(len(re.findall(r"^--mass-acc 20$", resolved, re.M)), 1, resolved)
            self.assertEqual(len(re.findall(r"^--mass-acc-ms1 7$", resolved, re.M)), 1, resolved)
            self.assertEqual(prov["search_mode"], "single_shot")
            ma = prov["result"]["mass_acc"]
            self.assertTrue(ma["measured"])
            self.assertEqual(ma["sop_floor"], {"--mass-acc": 20.0, "--mass-acc-ms1": 7.0})
            self.assertIn("max(measured, SOP)", ma["floor_note"])
            self.assertEqual(ma["value_file"], massacc)
            self.assertEqual(ma["evidence_file"], os.path.join(out, "mass_acc.json"))
            self.assertEqual(ma["documented"], {"--mass-acc-ms1": 7})

    def test_inline_single_shot_measures_before_the_search(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = self._run_search(d, *self._setup(d), "--allow-inline")
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(len(open(os.path.join(d, "argv.txt.probe")).read().splitlines()), 3)
            search = self._calls(d, "SEARCH")
            self.assertEqual(len(search), 1)
            self.assertIn("--mass-acc 20 --mass-acc-ms1 7", search[0])
            self.assertTrue(json.load(open(os.path.join(out, "search_provenance.json")))
                            ["result"]["ran"])

    def test_a_failed_measurement_stops_before_the_search(self):
        with tempfile.TemporaryDirectory() as d:
            raws = fx._cohort(d)
            for r in raws:                                  # no run reaches its MS2 line
                with open(r + ".log", "w") as fh:
                    fh.write(fx.no_ms2(fx.real_log(os.path.basename(r)[:-4])))
            p, out = self._run_search(d, *self._setup(d, runs=raws), "--allow-inline")
            self.assertNotEqual(p.returncode, 0)
            self.assertEqual(self._calls(d, "SEARCH"), [], "searched with nothing measured")
            self.assertIn("Ex01162023_10_TT33.raw", p.stdout + p.stderr)
            self.assertNotIn("Traceback", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(out, "params.resolved.cfg")))
            self.assertFalse(os.path.exists(os.path.join(out, "massacc.txt")))

    def test_a_pinned_cfg_is_not_probed(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = self._run_search(d, *self._setup(d, ms2_res=30000),
                                      "--sbatch", os.path.join(d, "job.sh"))
            self.assertEqual(p.returncode, 0, p.stderr)
            srch_sh = open(os.path.join(d, "job_2_search.sh")).read()
            self.assertNotIn("probe_window.py", srch_sh)
            self.assertNotIn("massacc.txt", srch_sh)

    def test_a_cfg_that_gets_xic_added_keeps_its_plan(self):
        """run_search.ensure_xic() searches a COPY of a cfg with no --xic. The copy must carry the
        rationale sidecar, or the plan is lost and the search silently goes back to auto mode."""
        with tempfile.TemporaryDirectory() as d:
            raws, fasta, cfg, tools, bundle = self._setup(d)
            kept = [ln for ln in open(cfg) if not ln.startswith(("--xic", "--mobilograms"))]
            with open(cfg, "w") as fh:
                fh.writelines(kept)
            p, out = self._run_search(d, raws, fasta, cfg, tools, bundle,
                                      "--sbatch", os.path.join(d, "job.sh"))
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("params_with_xic.cfg", p.stdout)
            self.assertIn("probe_window.py", open(os.path.join(d, "job_2_search.sh")).read())

    def test_a_cfg_without_a_sidecar_does_not_inherit_a_stale_one(self):
        """ensure_xic() copied the sidecar only when the source HAD one. A second search into the
        same --out with a hand-written cfg (no --xic, no mass accuracy, no sidecar) kept the first
        cfg's params_with_xic.cfg.rationale.json -- and its measure_with_diann plan, so the chain
        or the single-shot probe would measure mass accuracy and pin MS1 7 for a cfg that may not
        be an Orbitrap at all, instead of reporting the missing flags (mass_acc_unset)."""
        with tempfile.TemporaryDirectory() as d:
            _raws, _fasta, cfg, _tools, _bundle = self._setup(d)
            kept = [ln for ln in open(cfg) if not ln.startswith(("--xic", "--mobilograms"))]
            with open(cfg, "w") as fh:
                fh.writelines(kept)
            out = os.path.join(d, "out")
            aug = run_search.ensure_xic(cfg, out)
            self.assertIsNotNone(dp.mass_acc_measure_plan(aug))
            hand = os.path.join(d, "hand.cfg")
            with open(hand, "w") as fh:
                fh.writelines(kept)                       # the same flags, and no sidecar
            self.assertEqual(dp.parallel_safe(hand)["code"], "mass_acc_unset")
            self.assertEqual(run_search.ensure_xic(hand, out), aug)
            self.assertFalse(os.path.exists(aug + ".rationale.json"),
                             "the earlier cfg's rationale sidecar survived beside the new copy")
            self.assertIsNone(dp.mass_acc_measure_plan(aug))
            self.assertEqual(dp.parallel_safe(aug)["code"], "mass_acc_unset")

    def test_one_step_has_no_library_to_measure_against_and_says_so(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = self._run_search(d, *self._setup(d), "--one-step",
                                      "--sbatch", os.path.join(d, "job.sh"))
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertNotIn("probe_window.py", open(os.path.join(d, "job.sh")).read())
            self.assertIn("first run", p.stderr)
            ma = json.load(open(os.path.join(out, "search_provenance.json")))["result"]["mass_acc"]
            self.assertFalse(ma["measured"])
            self.assertIn("first run", ma["reason"])

    def test_the_search_job_gets_wall_clock_for_the_probe_it_runs(self):
        """The probe's --budget is diann_parallel.PROBE_BUDGET_S, sized against step 1b's OWN
        wall clock -- step 1b is a job of its own. Run inside the search job at the default 12 h,
        a probe that used its whole budget would have spent 3h50 of the search's wall before
        DIA-NN started, and SLURM would cut the search with the measurement done and nothing to
        show for it. The budget stays what step 1b proved; the job gets that much more wall."""
        base = run_search.SEARCH_WALL_HOURS

        def hours(script):
            (v,) = re.findall(r"^#SBATCH --time=(\d+):", open(script).read(), re.M)
            return int(v)

        with tempfile.TemporaryDirectory() as d:
            p, _ = self._run_search(d, *self._setup(d), "--sbatch", os.path.join(d, "job.sh"))
            self.assertEqual(p.returncode, 0, p.stderr)
            srch = os.path.join(d, "job_2_search.sh")
            self.assertIn("--budget %d" % dp.PROBE_BUDGET_S, open(srch).read())
            self.assertGreaterEqual(
                hours(srch) * 3600 - dp.PROBE_BUDGET_S, base * 3600,
                "the probe's budget is taken out of the search's own wall clock")
            self.assertEqual(hours(os.path.join(d, "job_1_lib.sh")), base,
                             "the library job does not probe and must not grow")
        with tempfile.TemporaryDirectory() as d:          # a pinned cfg: no probe, no extra wall
            args = self._setup(d, ms2_res=30000)
            p, _ = self._run_search(d, *args, "--sbatch", os.path.join(d, "job.sh"))
            self.assertEqual(p.returncode, 0, p.stderr)
            srch = os.path.join(d, "job_2_search.sh")
            self.assertNotIn("probe_window.py", open(srch).read())
            self.assertEqual(hours(srch), base)

    def _external_lib_cfg(self, d, cfg):
        """`cfg` rewritten to SEARCH an existing library instead of predicting one, keeping its
        measure_with_diann sidecar. This is a supported single-shot shape, not a broken cfg."""
        lib = os.path.join(d, "external.speclib")
        with open(lib, "w") as fh:
            fh.write("lib")
        out = os.path.join(d, "external.cfg")
        keep = [ln for ln in open(cfg).read().splitlines()
                if ln.split()[:1] not in (["--fasta-search"], ["--gen-spec-lib"],
                                          ["--predictor"])]
        with open(out, "w") as fh:
            fh.write("\n".join(keep + ["--lib " + lib]) + "\n")
        shutil.copy(cfg + ".rationale.json", out + ".rationale.json")
        return out

    def test_each_reason_there_is_nothing_to_measure_against_names_itself(self):
        """Three different cfg shapes reach the same "nothing was measured" note and they have
        three different fixes. It used to read "no library before it to measure against
        (--one-step, or a cfg that does not predict one) ... Drop --one-step to measure it" for
        all three -- wrong for a cfg that searches an external --lib, which HAS a library and has
        no --one-step to drop."""
        with tempfile.TemporaryDirectory() as d:
            raws, fasta, cfg, tools, bundle = self._setup(d)
            p, out = self._run_search(d, raws, fasta, cfg, tools, bundle, "--one-step",
                                      "--sbatch", os.path.join(d, "one.sh"))
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("--one-step was given", p.stderr)
            self.assertIn("Drop --one-step", p.stderr)

        with tempfile.TemporaryDirectory() as d:
            raws, fasta, cfg, tools, bundle = self._setup(d)
            ext = self._external_lib_cfg(d, cfg)
            self.assertIsNotNone(dp.mass_acc_measure_plan(ext), "the cfg still plans to measure")
            p, out = self._run_search(d, raws, fasta, ext, tools, bundle,
                                      "--sbatch", os.path.join(d, "ext.sh"))
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertNotIn("--one-step", p.stderr,
                             "an external-library cfg was told to drop a flag it never had")
            self.assertIn("--lib", p.stderr)
            self.assertIn("did not build", p.stderr)
            self.assertIn("--ms1-resolution", p.stderr)
            ma = json.load(open(os.path.join(out,
                                             "search_provenance.json")))["result"]["mass_acc"]
            self.assertFalse(ma["measured"])
            self.assertNotIn("--one-step", ma["reason"])

        with tempfile.TemporaryDirectory() as d:          # neither predicts nor supplies one
            raws, fasta, cfg, tools, bundle = self._setup(d)
            none_cfg = self._external_lib_cfg(d, cfg)
            keep = [ln for ln in open(none_cfg).read().splitlines()
                    if ln.split()[:1] != ["--lib"]]
            with open(none_cfg, "w") as fh:
                fh.write("\n".join(keep) + "\n")
            p, out = self._run_search(d, raws, fasta, none_cfg, tools, bundle,
                                      "--sbatch", os.path.join(d, "none.sh"))
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("neither predicts a library", p.stderr)
            self.assertNotIn("--one-step", p.stderr)
    def test_an_implausible_measurement_stops_the_search_job_and_searches_nothing(self):
        """The end of the path the band exists for: DIA-NN reports a tolerance, the probe refuses
        to pin it, and the search job stops there. Nothing is searched at 60 ppm, no massacc.txt
        or params.resolved.cfg is left for a resubmission to pick up, and the per-run evidence
        survives for whoever has to work out why (a wrong FASTA, a wrong species, a batch out of
        calibration -- none of which the median over three runs can see)."""
        with tempfile.TemporaryDirectory() as d:
            raws = fx._cohort(d)
            for r in raws:
                with open(r + ".log", "w") as fh:
                    fh.write(fx.ms2_as(fx.real_log(os.path.basename(r)[:-4]), 60))
            args = list(self._setup(d, runs=raws))
            p, out = self._run_search(d, *args, "--sbatch", os.path.join(d, "job.sh"))
            self.assertEqual(p.returncode, 0, p.stderr)
            env = self._env(d)
            r1 = subprocess.run(["bash", os.path.join(d, "job_1_lib.sh")], capture_output=True,
                                text=True, env=env, timeout=240)
            self.assertEqual(r1.returncode, 0, r1.stderr)
            r2 = subprocess.run(["bash", os.path.join(d, "job_2_search.sh")], capture_output=True,
                                text=True, env=env, timeout=240)
            self.assertNotEqual(r2.returncode, 0, "the search ran on a refused mass accuracy")
            self.assertIn("Refusing to pin", r2.stdout + r2.stderr)
            self.assertIn("Nothing was searched", r2.stdout + r2.stderr)
            for leftover in ("massacc.txt", "params.resolved.cfg", "params.resolved.cfg.tmp"):
                self.assertFalse(os.path.exists(os.path.join(out, leftover)), leftover)
            self.assertEqual(self._calls(d, "SEARCH"), [], "DIA-NN searched at 60 ppm")
            ev = json.load(open(os.path.join(out, "mass_acc.json")))
            self.assertIsNone(ev["mass_acc"])
            self.assertEqual([x["ms2_ppm"] for x in ev["probes"]], [60.0, 60.0, 60.0])
            self.assertTrue(any("outside 3-30 ppm" in r for r in ev["mass_acc_refused"]))

if __name__ == "__main__":
    unittest.main(verbosity=2)
