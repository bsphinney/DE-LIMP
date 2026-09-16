#!/usr/bin/env python3
"""
search_provenance.json must say which engine build ran, and say so when that is not the build
the workflow manifest asked for.

Two files name the engine version and nothing compared them:

  * workflow.manifest.json -- resolve_defaults.py's PIN (ENGINE_VERSIONS: diann 2.6.1). A request.
  * tools.json             -- what acquire_tools.sh resolved, next to the command run_search.py
                              executes verbatim. The build.

The FRAN pilot on 2026-09-16 ran DIA-NN 2.7.0 (banner "DIA-NN 2.7.0 Academia" on a compute node)
against a manifest pinning 2.6.1. run_search.py took `tools.versions or manifest.version`, so the
record happened to be right -- but it carried no trace of the disagreement, and when tools.json
had no entry the SAME expression quietly recorded the manifest's pin as if it had run (for any
engine: a Sage search under a DIA-NN manifest would have been stamped with DIA-NN's version).

tools.json is authoritative for `version`, because it describes the command that actually runs.
The manifest's pin is still recorded, both are compared, and a mismatch is printed and flagged.
"""
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import run_search  # noqa: E402


def record(engine, tools_versions, manifest_engine):
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        rec = run_search.engine_version_record(
            engine, {"versions": tools_versions, engine: "/opt/engine/bin"},
            {"engine": manifest_engine} if manifest_engine is not None else {})
    return rec, err.getvalue()


class EngineVersionRecordTests(unittest.TestCase):
    def test_pilot_mismatch_is_recorded_both_ways_and_warned(self):
        rec, err = record("diann", {"diann": "2.7.0"}, {"name": "diann", "version": "2.6.1"})
        self.assertEqual(rec["version"], "2.7.0")              # the build that runs
        self.assertEqual(rec["version_source"], "tools.json")
        self.assertEqual(rec["tools_engine_version"], "2.7.0")
        self.assertEqual(rec["manifest_engine_version"], "2.6.1")
        self.assertIs(rec["engine_version_mismatch"], True)
        self.assertIn("WARNING", err)
        self.assertIn("2.6.1", err)
        self.assertIn("2.7.0", err)

    def test_agreement_is_quiet(self):
        rec, err = record("diann", {"diann": "2.6.1"}, {"name": "diann", "version": "2.6.1"})
        self.assertEqual(rec["version"], "2.6.1")
        self.assertIs(rec["engine_version_mismatch"], False)
        self.assertEqual(err, "")

    def test_a_leading_v_is_not_a_mismatch(self):
        """Sage release tags are `v0.14.7`; resolve_defaults.py pins `0.14.7`."""
        rec, err = record("sage", {"sage": "v0.14.7"}, {"name": "sage", "version": "0.14.7"})
        self.assertIs(rec["engine_version_mismatch"], False)
        self.assertEqual(err, "")

    def test_a_placeholder_in_tools_json_falls_back_to_the_pin_LABELLED(self):
        """A tools.json written before acquire_tools.sh resolved `latest` says "latest" (and
        Sage from the conda env says "env"). Neither names a build. The pin is the best
        available answer, but the record must say it is a request nobody confirmed."""
        for placeholder in ("latest", "env", ""):
            rec, err = record("diann", {"diann": placeholder},
                              {"name": "diann", "version": "2.6.1"})
            self.assertEqual(rec["version"], "2.6.1", placeholder)
            self.assertNotEqual(rec["version_source"], "tools.json", placeholder)
            self.assertIn("manifest", rec["version_source"], placeholder)
            self.assertIs(rec["engine_version_mismatch"], False, placeholder)
            self.assertEqual(rec["tools_engine_version"], placeholder or None)
            self.assertIn("tools.json", err, placeholder)

    def test_another_engines_pin_is_never_borrowed(self):
        """The manifest pins ONE engine. `--engine sage` under a DIA-NN manifest used to be
        stamped with DIA-NN's version when tools.json had no Sage entry."""
        rec, _ = record("sage", {}, {"name": "diann", "version": "2.6.1"})
        self.assertIsNone(rec["version"])
        self.assertIsNone(rec["manifest_engine_version"])
        self.assertIs(rec["engine_version_mismatch"], False)

    def test_no_version_anywhere_is_null_and_said(self):
        rec, err = record("fragpipe", {}, {"name": "fragpipe"})
        self.assertIsNone(rec["version"])
        self.assertIsNone(rec["version_source"])
        self.assertIn("fragpipe", err)


class SearchProvenanceFileTests(unittest.TestCase):
    """End to end through main(): the fields land in search_provenance.json."""

    def _run(self, d, tools_diann_version, manifest_version):
        raw = os.path.join(d, "a.d")                   # a timsTOF-shaped input: no .NET needed
        os.makedirs(raw)
        cfg = os.path.join(d, "diann.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n--xic 10 --mobilograms\n")
        fasta = os.path.join(d, "db.fasta")
        with open(fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDEK\n")
        tools = os.path.join(d, "tools.json")
        with open(tools, "w") as fh:
            json.dump({"diann": "/opt/diann/diann-linux",
                       "versions": {"diann": tools_diann_version}}, fh)
        bundle = os.path.join(d, "workflow.manifest.json")
        with open(bundle, "w") as fh:
            json.dump({"acquisition": "DIA",
                       "engine": {"name": "diann", "version": manifest_version}}, fh)
        out = os.path.join(d, "search_out")
        # --sbatch only EMITS a script; nothing is executed. PATH without sbatch keeps the
        # routing single-shot and the login-node guard out of the way on any machine.
        env = dict(os.environ, PATH="/usr/bin:/bin")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                            "--tools", tools, "--bundle", bundle, "--params", cfg,
                            "--fasta", fasta, "--out", out, "--files", raw,
                            "--engine", "diann", "--sbatch", os.path.join(d, "job.sh")],
                           capture_output=True, text=True, env=env, timeout=120)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(out, "search_provenance.json")) as fh:
            return json.load(fh), r.stderr

    def test_the_pilot_mismatch_reaches_the_provenance_file(self):
        with tempfile.TemporaryDirectory() as d:
            prov, err = self._run(d, "2.7.0", "2.6.1")
            self.assertEqual(prov["version"], "2.7.0")
            self.assertEqual(prov["tools_engine_version"], "2.7.0")
            self.assertEqual(prov["manifest_engine_version"], "2.6.1")
            self.assertIs(prov["engine_version_mismatch"], True)
            self.assertEqual(prov["version_source"], "tools.json")
            self.assertIn("WARNING", err)

    def test_matching_versions_record_no_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            prov, _ = self._run(d, "2.6.1", "2.6.1")
            self.assertIs(prov["engine_version_mismatch"], False)
            self.assertEqual(prov["version"], "2.6.1")


if __name__ == "__main__":
    unittest.main()
