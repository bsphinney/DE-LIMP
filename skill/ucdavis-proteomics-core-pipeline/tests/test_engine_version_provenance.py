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
had no entry the SAME expression quietly recorded the manifest's pin as if it had run -- even
another engine's: acquire_tools.sh writes `versions` only for diann, sage and radiant, so a
FragPipe or AlphaDIA search under a DIA-NN manifest was stamped with DIA-NN's version.

`version` is only ever a build something CONFIRMS ran: tools.json `versions`, or failing that the
version the command itself names (a build folder like diann-2.6.0/, a .sif name, an image tag).
The manifest's pin is a request; it is recorded beside `version` and compared with it, but it is
never promoted to `version`, because fran_deposit.detect_engine() hands `version` -- and only
`version` -- to FRAN as the engine version. When nothing names the build, `version` is null and
`engine_version_mismatch` is null (not comparable), not false.
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


def record(engine, tools_versions, manifest_engine, cmd="/opt/engine/bin", **tools_extra):
    err = io.StringIO()
    tools = {"versions": tools_versions, engine: cmd, **tools_extra}
    with contextlib.redirect_stderr(err):
        rec = run_search.engine_version_record(
            engine, tools, {"engine": manifest_engine} if manifest_engine is not None else {})
    return rec, err.getvalue()


DIANN_261 = {"name": "diann", "version": "2.6.1"}


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

    def test_a_placeholder_never_promotes_the_pin_to_version(self):
        """A tools.json written before acquire_tools.sh resolved `latest` says "latest"; Sage
        from a conda env says "env". Neither names a build, and the command here names none
        either. The pin is what was ASKED for: recording it as `version` would send FRAN a
        version nothing checked (an old tools.json pointing at build_260/diann-2.6.0 was
        recorded as 2.6.1 this way). So: null, and not comparable."""
        for placeholder in ("latest", "env", ""):
            rec, err = record("diann", {"diann": placeholder}, DIANN_261)
            self.assertIsNone(rec["version"], placeholder)
            self.assertIsNone(rec["version_source"], placeholder)
            self.assertIsNone(rec["engine_version_mismatch"], placeholder)
            self.assertEqual(rec["manifest_engine_version"], "2.6.1", placeholder)
            self.assertEqual(rec["tools_engine_version"], placeholder or None)
            self.assertIn("null", err, placeholder)

    def test_an_old_latest_tools_json_takes_the_build_its_command_path_names(self):
        """An old acquire_tools.sh wrote "latest" beside a HIVE build path that names 2.6.0,
        under a manifest pinning 2.6.1. 2.6.0 is what runs."""
        rec, err = record("diann", {"diann": "latest"}, DIANN_261,
                          cmd="/quobyte/proteomics-grp/dia-nn/build_260/diann-2.6.0/diann-linux")
        self.assertEqual(rec["version"], "2.6.0")
        self.assertEqual(rec["version_source"], "command")
        self.assertEqual(rec["command_engine_versions"], ["2.6.0"])
        self.assertIs(rec["engine_version_mismatch"], True)
        self.assertIn("WARNING", err)

    def test_a_docker_image_tag_names_the_build(self):
        rec, _ = record("diann", {"diann": ""}, DIANN_261,
                        cmd="docker run --rm -v $PWD:/data proteomics-pipeline/diann:2.7.0 diann-linux")
        self.assertEqual(rec["version"], "2.7.0")
        self.assertIs(rec["engine_version_mismatch"], True)

    def test_a_hive_sif_names_the_build(self):
        rec, _ = record("diann", {"diann": "latest"}, DIANN_261,
                        cmd="apptainer exec --bind /quobyte:/quobyte "
                            "/quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif /diann-*/diann-linux")
        self.assertEqual(rec["version"], "2.3.0")

    def test_conda_env_sage_records_null_and_does_not_say_rerun_acquire_tools(self):
        """acquire_tools.sh writes "env" for EVERY sage on PATH, so re-running it cannot
        change the record -- telling the user to do that is advice that cannot work."""
        rec, err = record("sage", {"sage": "env"}, {"name": "sage", "version": "0.14.7"},
                          cmd="/home/u/.proteomics-pipeline/micromamba/envs/proteomics-pipeline/bin/sage")
        self.assertIsNone(rec["version"])
        self.assertIsNone(rec["engine_version_mismatch"])
        self.assertNotIn("acquire_tools", err)

    def test_an_unpinned_radiant_image_records_null(self):
        """tools.json `radiant` is only the runtime prefix; the image is `radiant_image`."""
        rec, _ = record("radiant", {"radiant": "latest"}, {"name": "radiant", "version": "2.3.3"},
                        cmd="docker run --rm", radiant_image="seerbio/radiant-fulcrum:latest")
        self.assertIsNone(rec["version"])
        rec, _ = record("radiant", {"radiant": ""}, {"name": "radiant", "version": "2.3.3"},
                        cmd="docker run --rm", radiant_image="seerbio/radiant-fulcrum:2.3.3")
        self.assertEqual(rec["version"], "2.3.3")
        self.assertIs(rec["engine_version_mismatch"], False)

    def test_tools_json_contradicting_its_own_command_records_null(self):
        """tools.json says 2.6.1 but the path it runs is diann-2.6.0/. One of them is wrong
        and nothing here can tell which: say so, record neither as `version`."""
        rec, err = record("diann", {"diann": "2.6.1"}, DIANN_261,
                          cmd="/opt/dia-nn/build_260/diann-2.6.0/diann-linux")
        self.assertIsNone(rec["version"])
        self.assertIsNone(rec["engine_version_mismatch"])
        self.assertEqual(rec["tools_engine_version"], "2.6.1")
        self.assertEqual(rec["command_engine_versions"], ["2.6.0"])
        self.assertIn("WARNING", err)

    def test_tools_json_agreeing_with_its_command_is_quiet(self):
        rec, err = record("diann", {"diann": "2.6.1"}, DIANN_261,
                          cmd="/quobyte/proteomics-grp/dia-nn/build_261/diann-2.6.1/diann-linux")
        self.assertEqual(rec["version"], "2.6.1")
        self.assertEqual(rec["version_source"], "tools.json")
        self.assertIs(rec["engine_version_mismatch"], False)
        self.assertEqual(err, "")

    def test_the_request_keyed_sage_cache_dir_is_not_read_as_a_version(self):
        """acquire_sage caches under sage/<REQUEST>/: sage/0.14.6/sage can hold another
        release (acquire_tools.sh notes that case). Only a release folder name counts."""
        rec, _ = record("sage", {"sage": ""}, {"name": "sage", "version": "0.14.7"},
                        cmd="/home/u/.proteomics-pipeline/tools/sage/0.14.6/sage")
        self.assertIsNone(rec["version"])
        self.assertEqual(rec["command_engine_versions"], [])

    def test_mismatch_fix_it_names_this_tools_root_and_platform(self):
        """The pilot ran with --tools .../pilot_tools/tools.json. A fix-it that omits the root
        rewrites ~/.proteomics-pipeline/tools/tools.json instead, and the WARNING returns."""
        rec, err = record("diann", {"diann": "2.7.0"}, DIANN_261,
                          cmd="/quobyte/proteomics-grp/fran/engines/pilot_tools/diann/2.7.0/diann-2.7.0/diann-linux",
                          tools_root="/quobyte/proteomics-grp/fran/engines/pilot_tools",
                          platform_class="hpc")
        self.assertIs(rec["engine_version_mismatch"], True)
        self.assertIn("PIN_ENGINE=diann PIN_VERSION=2.6.1 bash scripts/acquire_tools.sh hpc "
                      "/quobyte/proteomics-grp/fran/engines/pilot_tools", err)

    def test_mac_docker_mismatch_fix_it_rebuilds_the_image(self):
        """On macOS acquire_tools.sh only wraps $DIANN_DOCKER_IMAGE, so re-running it cannot
        change the version; the image has to be built for the pin."""
        rec, err = record("diann", {"diann": "2.7.0"}, DIANN_261,
                          cmd="docker run --rm -v $PWD:/data proteomics-pipeline/diann:2.7.0 diann-linux",
                          platform_class="mac", tools_root="/Users/u/.proteomics-pipeline/tools")
        self.assertIs(rec["engine_version_mismatch"], True)
        self.assertIn("build_diann_docker.sh 2.6.1", err)

    def test_another_engines_pin_is_never_borrowed(self):
        """The manifest pins ONE engine. tools.json has no `versions` entry for AlphaDIA or
        FragPipe, and `--engine alphadia` under a DIA-NN manifest used to be stamped with
        DIA-NN's version."""
        rec, _ = record("alphadia", {"diann": "2.6.1", "sage": "latest"},
                        {"name": "diann", "version": "2.6.1"},
                        cmd="apptainer exec --bind /quobyte:/quobyte "
                            "/quobyte/proteomics-grp/apptainers/alphadia.sif alphadia")
        self.assertIsNone(rec["version"])
        self.assertIsNone(rec["manifest_engine_version"])
        self.assertIsNone(rec["engine_version_mismatch"])

    def test_no_version_anywhere_is_null_and_said(self):
        rec, err = record("fragpipe", {}, {"name": "fragpipe"})
        self.assertIsNone(rec["version"])
        self.assertIsNone(rec["version_source"])
        self.assertIn("fragpipe", err)


class SearchProvenanceFileTests(unittest.TestCase):
    """End to end through main(): the fields land in search_provenance.json."""

    def _run(self, d, tools_diann_version, manifest_version, diann_cmd="/opt/diann/diann-linux"):
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
            json.dump({"diann": diann_cmd,
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

    def test_what_fran_deposit_forwards_is_never_the_unconfirmed_pin(self):
        """fran_deposit.detect_engine() passes `version` -- and nothing else from this record --
        to FRAN as engine_version. Check at that consumer, not only in the file."""
        import fran_deposit
        with tempfile.TemporaryDirectory() as d:
            self._run(d, "latest", "2.6.1")                     # names no build anywhere
            eng, ver, _src = fran_deposit.detect_engine(os.path.join(d, "search_out"))
            self.assertEqual(eng, "diann")
            self.assertIsNone(ver)
        with tempfile.TemporaryDirectory() as d:
            self._run(d, "latest", "2.6.1",
                      diann_cmd="/quobyte/proteomics-grp/dia-nn/build_260/diann-2.6.0/diann-linux")
            _eng, ver, _src = fran_deposit.detect_engine(os.path.join(d, "search_out"))
            self.assertEqual(ver, "2.6.0")


if __name__ == "__main__":
    unittest.main()
