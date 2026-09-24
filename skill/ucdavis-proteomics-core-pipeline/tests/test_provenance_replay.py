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


class OrbitrapResolutionReplayTests(unittest.TestCase):
    """reproduce.sh re-derived the defaults WITHOUT the Orbitrap resolution, so replaying a
    60k/15k Fusion Lumos run reclassified it as orbitrap_generic and derived different mass
    accuracy (review 2026-09-24). The replay now passes the recorded resolution AND its source,
    and the rebuilt manifest must equal the original."""

    def test_replay_carries_the_resolution_and_its_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = write_env(tmp)
            wf = os.path.join(tmp, "wf")
            base = ["--acquisition", "DIA", "--instrument", "Orbitrap Fusion Lumos",
                    "--engine", "diann", "--organism-taxid", "9606", "--env", env]
            subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                            *base, "--ms1-resolution", "60000", "--ms2-resolution", "15000",
                            "--resolution-source", "detected", "--dest", wf],
                           capture_output=True, text=True, check=True)
            out = os.path.join(tmp, "repro")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                                "--outdir", out, "--workflow-manifest",
                                os.path.join(wf, "workflow.manifest.json"),
                                "--engine", "diann", "--acquisition", "DIA",
                                "--instrument", "Orbitrap Fusion Lumos",
                                "--organism-taxid", "9606"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(os.path.join(out, "reproduce.sh")) as fh:
                sh = fh.read()
            self.assertIn("--ms1-resolution 60000 --ms2-resolution 15000 "
                          "--resolution-source detected", sh)
            # replay exactly those flags
            replay = os.path.join(tmp, "replay")
            subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                            *base, "--ms1-resolution", "60000", "--ms2-resolution", "15000",
                            "--resolution-source", "detected", "--dest", replay],
                           capture_output=True, text=True, check=True)

            def load(p):
                with open(p) as fh:
                    d = json.load(fh)
                d.pop("output", None)
                (d.get("search") or {}).pop("params_file", None)
                return d
            a = load(os.path.join(wf, "workflow.manifest.json"))
            b = load(os.path.join(replay, "workflow.manifest.json"))
            self.assertEqual((a.get("resolution") or {}).get("source"), "detected")
            self.assertEqual([k for k in sorted(set(a) | set(b)) if a.get(k) != b.get(k)], [])


class IonTrapReplayTests(unittest.TestCase):
    """An Orbitrap-MS1 / ion-trap-MS2 run has no MS2 resolution; --ms2-analyzer ITMS is what
    selects its tolerances, so the replay must carry it."""

    def test_replay_carries_the_ion_trap_analyzer(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = write_env(tmp)
            wf = os.path.join(tmp, "wf")
            subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                            "--acquisition", "DDA", "--instrument", "Orbitrap Fusion Lumos",
                            "--engine", "sage", "--organism-taxid", "9606", "--env", env,
                            "--ms1-resolution", "120000", "--ms2-analyzer", "ITMS",
                            "--resolution-source", "detected", "--dest", wf],
                           capture_output=True, text=True, check=True)
            out = os.path.join(tmp, "repro")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                                "--outdir", out, "--workflow-manifest",
                                os.path.join(wf, "workflow.manifest.json"),
                                "--engine", "sage", "--acquisition", "DDA",
                                "--instrument", "Orbitrap Fusion Lumos",
                                "--organism-taxid", "9606"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(os.path.join(out, "reproduce.sh")) as fh:
                sh = fh.read()
            self.assertIn("--ms1-resolution 120000 --ms2-analyzer ITMS "
                          "--resolution-source detected", sh)
            self.assertNotIn("--ms2-resolution", sh)


class DigestionEnzymeReplayTests(unittest.TestCase):
    """fetch_fasta.py --enzyme decides which protease contaminant entries stay Cont_ when they
    match a target protein; reproduce.sh must rebuild the database with the same enzymes, and
    fall back to the default (trypsin,lysc) for sidecars written before --enzyme existed."""

    def _repro(self, tmp, fasta_info):
        out, wf = build_bundle(tmp)
        out2 = os.path.join(tmp, "repro2")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                            "--outdir", out2, "--workflow-manifest",
                            os.path.join(wf, "workflow.manifest.json"),
                            "--engine", "diann", "--acquisition", "DIA",
                            "--instrument", "timsTOF HT", "--organism-taxid", "10090",
                            "--fasta-info", json.dumps(fasta_info)],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(out2, "reproduce.sh")) as fh:
            return fh.read()

    def test_non_default_enzyme_is_replayed(self):
        with tempfile.TemporaryDirectory() as tmp:
            sh = self._repro(tmp, {"proteome": "UP000000589", "content_used": "one_per_gene",
                                   "contaminant_set": "universal",
                                   "digestion_enzymes_used": ["gluc"]})
            self.assertIn("--enzyme gluc", sh)

    def test_old_sidecar_replays_the_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            sh = self._repro(tmp, {"proteome": "UP000000589", "content_used": "one_per_gene",
                                   "contaminant_set": "universal"})
            self.assertIn("--enzyme trypsin,lysc", sh)


class SageSelfReportedVersionTests(unittest.TestCase):
    """environment/versions.txt must not present `sage --version` as THE Sage version.

    The v0.14.7 Sage release binary prints "sage 0.14.6" (measured on HIVE 2026-09-16), which
    is why acquire_tools.sh reads the release from the tarball instead. This file wrote the
    binary's self-report under the flat name `sage_version`, so one run produced two records
    that disagreed -- versions.txt saying 0.14.6 and search_provenance.json saying 0.14.7 --
    with nothing on either to say which was the release and which was the binary talking.
    """

    def _versions(self, tools=None):
        with tempfile.TemporaryDirectory() as tmp:
            sage = os.path.join(tmp, "sage")
            with open(sage, "w") as fh:
                fh.write("#!/bin/sh\necho 'sage 0.14.6'\n")
            os.chmod(sage, 0o755)
            setup_json = os.path.join(tmp, "setup.json")
            with open(setup_json, "w") as fh:
                json.dump({"sage": sage}, fh)
            argv = [sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                    "--outdir", os.path.join(tmp, "repro"), "--engine", "sage",
                    "--setup-json", setup_json]
            if tools is not None:
                tj = os.path.join(tmp, "tools.json")
                with open(tj, "w") as fh:
                    json.dump(tools, fh)
                argv += ["--tools-json", tj]
            r = subprocess.run(argv, capture_output=True, text=True, timeout=300)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            with open(os.path.join(tmp, "repro", "environment", "versions.txt")) as fh:
                return json.load(fh)

    def test_the_binarys_self_report_is_not_recorded_as_the_sage_version(self):
        v = self._versions()
        self.assertNotIn("sage_version", v, "an unqualified name reads as the release")
        self.assertIn("sage 0.14.6", v["sage_self_reported_version"]["value"])

    def test_it_carries_the_caveat_that_makes_the_disagreement_readable(self):
        v = self._versions({"versions": {"sage": "0.14.7"}, "sage": "/opt/sage/sage"})
        self.assertEqual(v["tools_versions"]["sage"], "0.14.7")     # the release of record
        note = v["sage_self_reported_version"]["note"]
        self.assertIn("0.14.6", note)
        self.assertIn("tools_versions.sage", note)


if __name__ == "__main__":
    unittest.main(verbosity=2)
