#!/usr/bin/env python3
"""setup.json's next_steps: the fix for each item that is not ready, as one sentence for this
platform, so the agent relays it instead of deriving it from ready_for (SKILL.md step 0).

Keys de, limpa, dia, dda, thermo_raw, resolution_reader; a ready item has no entry, so a ready
machine has {}. Each sentence says when it matters (DIA data only, .raw input only) and reuses the
fix setup.sh already gives elsewhere (limpa's LIMPA_FIX, nanoparquet.fix, diann.note,
thermo_raw_reader.note) rather than a second wording of it; each is also printed to stderr as a
"[setup] next:" line. Built through test_setup_limpa's stand-ins (stub Rscript + micromamba, no
conda, no network), with `setup.sh --check`.
"""
import json
import os
import platform
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import test_setup_limpa as limpa_t                    # noqa: E402
import test_thermo_raw_detection as trfp              # noqa: E402
from test_trfp_dotnet_discovery import _exe           # noqa: E402

IS_MAC = platform.system() == "Darwin"
FAKES = os.path.join(HERE, "fixtures", "fake_rawfilereader")
RAWFILEREADER_DLLS = ("ThermoFisher.CommonCore.RawFileReader.dll",
                      "ThermoFisher.CommonCore.Data.dll")
KEYS = {"de", "limpa", "dia", "dda", "thermo_raw", "resolution_reader"}


def folded(text):
    return " ".join(text.split()).rstrip(".")


class SetupNextSteps(limpa_t.SetupLimpaGate):
    def ready_env(self, limpa="1.4.1", parser=True):
        """An env where everything is ready: limpa, nanoparquet, Sage, a parser with the
        RawFileReader DLLs beside it (resolution read through the fake pythonnet), and docker
        on PATH (DIA-NN's only route on macOS)."""
        self.prebuilt_env(limpa)
        env_bin = os.path.join(self.prefix, "bin")
        _exe(os.path.join(env_bin, "sage"), "#!/bin/sh\nexit 0\n")
        if parser:
            _exe(os.path.join(env_bin, "ThermoRawFileParser"),
                 f'#!/bin/sh\nexec "{sys.executable}" "{trfp.FAKE}" "$@"\n')
            for name in RAWFILEREADER_DLLS:
                open(os.path.join(env_bin, name), "w").close()
        _exe(os.path.join(self.sys, "docker"), "#!/bin/sh\nexit 0\n")
        py_s = _exe(os.path.join(self.d, "py_S"), f'#!/bin/sh\nexec "{sys.executable}" -S "$@"\n')
        return {"THERMO_RESOLUTION_PYTHON": py_s, "PYTHONPATH": FAKES}

    def check(self, **env):
        r, s = self.run_setup("--check", **env)
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertLessEqual(set(s["next_steps"]), KEYS)
        # the keys other scripts read are still there, beside it
        self.assertEqual(set(s["ready_for"]), {"de", "dia", "dda", "thermo_raw"})
        for key in ("limpa", "de_reader", "thermo_raw_reader", "notes", "set_annotation"):
            self.assertIn(key, s)
        # one stderr line per step, the same sentence
        lines = [ln for ln in r.stderr.splitlines() if ln.startswith("[setup] next:")]
        self.assertEqual(lines, [f"[setup] next: {k}: {v}" for k, v in s["next_steps"].items()])
        return r, s

    def test_all_ready_is_empty(self):
        r, s = self.check(**self.ready_env())
        self.assertEqual(s["ready_for"], {"de": True, "dia": True, "dda": True, "thermo_raw": True})
        self.assertIs(s["limpa"]["ok"], True)
        self.assertIs(s["thermo_raw_reader"]["resolution_reader"]["ready"], True,
                      s["thermo_raw_reader"]["resolution_reader"])
        self.assertEqual(s["next_steps"], {})
        self.assertNotIn("[setup] next:", r.stderr)

    @unittest.skipUnless(IS_MAC, "the macOS branch: DIA-NN only through Docker")
    def test_dia_on_a_mac_without_docker_is_build_diann_docker(self):
        env = self.ready_env()
        os.remove(os.path.join(self.sys, "docker"))
        _r, s = self.check(**env)
        self.assertIs(s["ready_for"]["dia"], False)
        self.assertEqual(set(s["next_steps"]), {"dia"})
        step = s["next_steps"]["dia"]
        self.assertTrue(step.startswith("Only for DIA data"), step)
        self.assertIn(f"bash {limpa_t.setup_t.SCRIPTS}/build_diann_docker.sh", step)
        self.assertIn("never by falling back to a DDA engine", step)
        self.assertIn(s["diann"]["note"], step)                  # its fix, not a second wording

    def test_thermo_raw_relays_the_readers_own_note(self):
        _r, s = self.check(**self.ready_env(parser=False))
        self.assertIs(s["ready_for"]["thermo_raw"], False)
        step = s["next_steps"]["thermo_raw"]
        self.assertTrue(step.startswith("Only if the input has Thermo .raw files"), step)
        self.assertIn(folded(s["thermo_raw_reader"]["note"]), step)
        self.assertIn("setup.sh --check until ready_for.thermo_raw is true", step)
        # no parser, so no DLLs beside it: the resolution cannot be read either, and its note says so
        rr = s["thermo_raw_reader"]["resolution_reader"]
        self.assertIs(rr["ready"], False)
        self.assertIn(folded(rr["note"]), s["next_steps"]["resolution_reader"])
        self.assertTrue(s["next_steps"]["resolution_reader"].startswith(
            "Only for Thermo Orbitrap .raw input"))
        self.assertEqual(set(s["next_steps"]), {"thermo_raw", "resolution_reader"})

    def test_an_old_limpa_says_ask_the_director_first(self):
        r, s = self.check(**self.ready_env(limpa="1.2.5"))
        self.assertIs(s["limpa"]["ok"], False)
        self.assertIs(s["ready_for"]["de"], True)                 # limpa 1.2.x reads with arrow
        self.assertEqual(set(s["next_steps"]), {"limpa"})
        step = s["next_steps"]["limpa"]
        self.assertTrue(step.startswith("limpa 1.2.5 is older than 1.4.0"), step)
        self.assertIn("Core staff ask the Core's director first", step)
        # LIMPA_FIX, the fix the ERROR line prints
        fix = r.stderr.split("[setup]        Fix: ", 1)[1].splitlines()[0]
        self.assertIn(fix, step)

    def test_no_limpa_at_all_is_an_install_not_an_upgrade(self):
        """Nothing has run on a limpa that is not there: no director to ask."""
        _r, s = self.check(**self.ready_env(limpa=""))
        self.assertEqual(set(s["next_steps"]), {"de", "limpa"})
        self.assertTrue(s["next_steps"]["limpa"].startswith("limpa is not installed"))
        self.assertNotIn("director", s["next_steps"]["limpa"])
        self.assertIn("--method maxlfq can still run", s["next_steps"]["de"])

    def test_de_without_nanoparquet_gives_its_fix_and_says_maxlfq_runs(self):
        _r, s = self.check(**dict(self.ready_env(), **self.np_missing()))
        self.assertIs(s["ready_for"]["de"], False)
        self.assertEqual(set(s["next_steps"]), {"de"})
        step = s["next_steps"]["de"]
        self.assertIn(s["nanoparquet"]["fix"], step)
        self.assertEqual(step.count("no DE can run"), 1, step)     # the note's own wording, once
        self.assertIn("--method maxlfq can still run", step)
        self.assertNotIn("login node", step)                      # not a cluster

    def test_on_a_cluster_the_fix_never_runs_on_the_login_node(self):
        env = dict(self.ready_env(), **self.np_missing(), SLURM_JOB_ID="1")
        _r, s = self.check(**env)
        step = s["next_steps"]["de"]
        self.assertIn("from a login node it runs itself on a compute node", step)
        self.assertIn("by hand inside a job, never on a login node: " + s["nanoparquet"]["fix"], step)

    def test_nothing_built_yet(self):
        _r, s = self.check()
        expected = {"de", "limpa", "dda", "thermo_raw"} | ({"dia"} if IS_MAC else set())
        self.assertEqual(set(s["next_steps"]), expected)
        self.assertIn("is not built yet", s["next_steps"]["limpa"])
        self.assertNotIn("director", s["next_steps"]["limpa"])     # no env yet: nothing shared
        self.assertIn("Sage and R not found", s["next_steps"]["dda"])
        self.assertIn(folded(s["thermo_raw_reader"]["note"]), s["next_steps"]["thermo_raw"])
        if IS_MAC:
            self.assertIn("msconvert is Linux-only", s["next_steps"]["dda"])

    @unittest.skipUnless(os.access("/bin/bash", os.X_OK), "no /bin/bash")
    def test_the_system_bash_writes_the_same_json(self):
        """macOS's /bin/bash is 3.2: no associative arrays, and set -u with empty arrays."""
        env = self.ready_env(parser=False)
        self.bash = "/bin/bash"
        _r, s = self.check(**env)
        self.assertEqual(set(s["next_steps"]), {"thermo_raw", "resolution_reader"})
        with open(os.path.join(self.home, "setup.json")) as fh:
            json.load(fh)


# the parent's own tests run in test_setup_limpa.py; not again here
for _name in [n for n in dir(SetupNextSteps) if n.startswith("test_")]:
    if _name not in SetupNextSteps.__dict__:
        setattr(SetupNextSteps, _name, None)
del _name

if __name__ == "__main__":
    unittest.main()
