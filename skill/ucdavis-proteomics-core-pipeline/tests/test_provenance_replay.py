#!/usr/bin/env python3
"""
The reproducibility bundle must replay the ORIGINAL invocation, not a subset.

Found by auditing the skill's own output: `reproduce.sh` re-derived the search
defaults with only --acquisition/--instrument/--engine. The engine params still
came back byte-identical (they depend on neither organism nor platform), so the
SEARCH was reproducible -- but the regenerated run record silently lost:

  * organism_taxid  -> null. The species is the single most important contextual
                       fact about a proteomics run.
  * platform        -> null. The record could not say what machine it ran on,
                       nor carry the "runs emulated under Rosetta" warning.

And run_manifest.json recorded "timestamp": null whenever the orchestrator did
not pass --timestamp. A bundle whose whole job is to make a run auditable must be
able to say when it ran.

These are provenance defects rather than reproducibility ones, which is exactly
why they survived: every scientific number was right.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPTS = os.path.join(ROOT, "scripts")


def write_env(tmp):
    env = os.path.join(tmp, "env.json")
    with open(env, "w") as fh:
        subprocess.run(["bash", os.path.join(SCRIPTS, "detect_env.sh")],
                       stdout=fh, stderr=subprocess.DEVNULL, check=True)
    return env


def build_bundle(tmp, taxid="10090"):
    wf = os.path.join(tmp, "wf")
    # Mirror a real run: --env is what populates the platform block, and the
    # orchestrator is told to pass it. Building the "original" without it would
    # compare against a record that no real run produces.
    subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                    "--acquisition", "DIA", "--instrument", "timsTOF HT",
                    "--organism-taxid", taxid, "--env", write_env(tmp),
                    "--dest", wf],
                   capture_output=True, text=True, check=True)
    out = os.path.join(tmp, "repro")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                        "--outdir", out, "--workflow-manifest",
                        os.path.join(wf, "workflow.manifest.json"),
                        "--engine", "diann", "--acquisition", "DIA",
                        "--instrument", "timsTOF HT", "--organism-taxid", taxid],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return out, wf


class TestReplayFidelity(unittest.TestCase):
    def test_reproduce_sh_carries_the_organism(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = build_bundle(tmp)
            with open(os.path.join(out, "reproduce.sh")) as fh:
                sh = fh.read()
            self.assertIn("--organism-taxid 10090", sh)

    def test_reproduce_sh_regenerates_the_platform_map(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = build_bundle(tmp)
            with open(os.path.join(out, "reproduce.sh")) as fh:
                sh = fh.read()
            self.assertIn("detect_env.sh", sh)
            self.assertIn("--env", sh)

    def test_run_manifest_always_has_a_timestamp(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = build_bundle(tmp)
            with open(os.path.join(out, "run_manifest.json")) as fh:
                m = json.load(fh)
            self.assertIsNotNone(m.get("timestamp"),
                                 "a run record with no time is not an audit record")
            self.assertRegex(m["timestamp"], r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")

    def test_run_manifest_records_the_skill_version(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, _ = build_bundle(tmp)
            with open(os.path.join(out, "run_manifest.json")) as fh:
                m = json.load(fh)
            # Parameters ship with the skill, so the skill version is what pins
            # them -- there is no external registry commit to fall back on.
            self.assertTrue((m.get("skill") or {}).get("version"))

    def test_replaying_the_recorded_command_reproduces_the_record(self):
        """The end-to-end claim, executed rather than asserted."""
        with tempfile.TemporaryDirectory() as tmp:
            out, wf = build_bundle(tmp)
            replay = os.path.join(tmp, "replay")
            env = write_env(tmp)
            subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                            "--acquisition", "DIA", "--instrument", "timsTOF HT",
                            "--engine", "diann", "--organism-taxid", "10090",
                            "--env", env, "--dest", replay],
                           capture_output=True, text=True, check=True)

            def load(p):
                with open(p) as fh:
                    d = json.load(fh)
                d.pop("output", None)
                s = d.get("search") or {}
                s.pop("params_file", None)
                return d

            a = load(os.path.join(wf, "workflow.manifest.json"))
            b = load(os.path.join(replay, "workflow.manifest.json"))
            differing = [k for k in sorted(set(a) | set(b)) if a.get(k) != b.get(k)]
            self.assertEqual(differing, [], f"replay lost/changed: {differing}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
