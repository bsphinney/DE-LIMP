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
the mismatch is null (not comparable), not false.

The record has the shape of `scan_window`, its neighbour in search_provenance.json: ONE object,
`engine_version`, with a `value` and a `source` that says in words where the value came from
("unknown -- ..." when there is none), and every input kept beside them as written. Top-level
`version` stays, because fran_deposit.py forwards it to FRAN, and always equals
`engine_version.value`.
"""
import contextlib
import io
import json
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
        self.assertEqual(rec["value"], "2.7.0")                # the build that runs
        self.assertTrue(rec["source"].startswith("tools.json"), rec["source"])
        self.assertEqual(rec["tools_json"], "2.7.0")
        self.assertEqual(rec["manifest_pin"], "2.6.1")
        self.assertIs(rec["mismatch"], True)
        self.assertIn("WARNING", err)
        self.assertIn("2.6.1", err)
        self.assertIn("2.7.0", err)

    def test_agreement_is_quiet(self):
        rec, err = record("diann", {"diann": "2.6.1"}, {"name": "diann", "version": "2.6.1"})
        self.assertEqual(rec["value"], "2.6.1")
        self.assertIs(rec["mismatch"], False)
        self.assertEqual(err, "")

    def test_a_leading_v_is_not_a_mismatch(self):
        """Sage release tags are `v0.14.7`; resolve_defaults.py pins `0.14.7`."""
        rec, err = record("sage", {"sage": "v0.14.7"}, {"name": "sage", "version": "0.14.7"})
        self.assertIs(rec["mismatch"], False)
        self.assertEqual(err, "")

    def test_a_placeholder_never_promotes_the_pin_to_version(self):
        """A tools.json written before acquire_tools.sh resolved `latest` says "latest"; Sage
        from a conda env says "env". Neither names a build, and the command here names none
        either. The pin is what was ASKED for: recording it as `version` would send FRAN a
        version nothing checked (an old tools.json pointing at build_260/diann-2.6.0 was
        recorded as 2.6.1 this way). So: null, and not comparable."""
        for placeholder in ("latest", "env", ""):
            rec, err = record("diann", {"diann": placeholder}, DIANN_261)
            self.assertIsNone(rec["value"], placeholder)
            self.assertTrue(rec["source"].startswith("unknown -- "), rec["source"])
            self.assertNotIn("2.6.1", rec["source"], "the pin must not read as the source")
            self.assertIsNone(rec["mismatch"], placeholder)
            self.assertEqual(rec["manifest_pin"], "2.6.1", placeholder)
            self.assertEqual(rec["tools_json"], placeholder or None)
            self.assertIn("null", err, placeholder)

    def test_an_old_latest_tools_json_takes_the_build_its_command_path_names(self):
        """An old acquire_tools.sh wrote "latest" beside a HIVE build path that names 2.6.0,
        under a manifest pinning 2.6.1. 2.6.0 is what runs."""
        rec, err = record("diann", {"diann": "latest"}, DIANN_261,
                          cmd="/quobyte/proteomics-grp/dia-nn/build_260/diann-2.6.0/diann-linux")
        self.assertEqual(rec["value"], "2.6.0")
        self.assertTrue(rec["source"].startswith("named by the command"), rec["source"])
        self.assertIn("build_260/diann-2.6.0", rec["source"])
        self.assertEqual(rec["named_by_command"], ["2.6.0"])
        self.assertIs(rec["mismatch"], True)
        self.assertIn("WARNING", err)

    def test_a_docker_image_tag_names_the_build(self):
        rec, _ = record("diann", {"diann": ""}, DIANN_261,
                        cmd="docker run --rm -v $PWD:/data proteomics-pipeline/diann:2.7.0 diann-linux")
        self.assertEqual(rec["value"], "2.7.0")
        self.assertIs(rec["mismatch"], True)

    def test_a_hive_sif_names_the_build(self):
        rec, _ = record("diann", {"diann": "latest"}, DIANN_261,
                        cmd="apptainer exec --bind /quobyte:/quobyte "
                            "/quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif /diann-*/diann-linux")
        self.assertEqual(rec["value"], "2.3.0")

    def test_conda_env_sage_records_null_and_does_not_say_rerun_acquire_tools(self):
        """acquire_tools.sh writes "env" for EVERY sage on PATH, so re-running it cannot
        change the record -- telling the user to do that is advice that cannot work."""
        rec, err = record("sage", {"sage": "env"}, {"name": "sage", "version": "0.14.7"},
                          cmd="/home/u/.proteomics-pipeline/micromamba/envs/proteomics-pipeline/bin/sage")
        self.assertIsNone(rec["value"])
        self.assertIsNone(rec["mismatch"])
        self.assertNotIn("acquire_tools", err)

    def test_an_unpinned_radiant_image_records_null(self):
        """tools.json `radiant` is only the runtime prefix; the image is `radiant_image`."""
        rec, _ = record("radiant", {"radiant": "latest"}, {"name": "radiant", "version": "2.3.3"},
                        cmd="docker run --rm", radiant_image="seerbio/radiant-fulcrum:latest")
        self.assertIsNone(rec["value"])
        rec, _ = record("radiant", {"radiant": ""}, {"name": "radiant", "version": "2.3.3"},
                        cmd="docker run --rm", radiant_image="seerbio/radiant-fulcrum:2.3.3")
        self.assertEqual(rec["value"], "2.3.3")
        self.assertIs(rec["mismatch"], False)

    def test_tools_json_contradicting_its_own_command_records_null(self):
        """tools.json says 2.6.1 but the path it runs is diann-2.6.0/. One of them is wrong
        and nothing here can tell which: say so, record neither as `version`."""
        rec, err = record("diann", {"diann": "2.6.1"}, DIANN_261,
                          cmd="/opt/dia-nn/build_260/diann-2.6.0/diann-linux")
        self.assertIsNone(rec["value"])
        self.assertIsNone(rec["mismatch"])
        self.assertEqual(rec["tools_json"], "2.6.1")
        self.assertEqual(rec["named_by_command"], ["2.6.0"])
        self.assertTrue(rec["source"].startswith("unknown -- "), rec["source"])
        self.assertIn("2.6.0", rec["source"])
        self.assertIn("WARNING", err)

    def test_tools_json_agreeing_with_its_command_is_quiet(self):
        rec, err = record("diann", {"diann": "2.6.1"}, DIANN_261,
                          cmd="/quobyte/proteomics-grp/dia-nn/build_261/diann-2.6.1/diann-linux")
        self.assertEqual(rec["value"], "2.6.1")
        self.assertTrue(rec["source"].startswith("tools.json"), rec["source"])
        self.assertIs(rec["mismatch"], False)
        self.assertEqual(err, "")

    def test_the_request_keyed_sage_cache_dir_is_not_read_as_a_version(self):
        """acquire_sage caches under sage/<REQUEST>/: sage/0.14.6/sage can hold another
        release (acquire_tools.sh notes that case). Only a release folder name counts."""
        rec, _ = record("sage", {"sage": ""}, {"name": "sage", "version": "0.14.7"},
                        cmd="/home/u/.proteomics-pipeline/tools/sage/0.14.6/sage")
        self.assertIsNone(rec["value"])
        self.assertEqual(rec["named_by_command"], [])

    def test_mismatch_fix_it_names_this_tools_root_and_platform(self):
        """The pilot ran with --tools .../pilot_tools/tools.json. A fix-it that omits the root
        rewrites ~/.proteomics-pipeline/tools/tools.json instead, and the WARNING returns."""
        rec, err = record("diann", {"diann": "2.7.0"}, DIANN_261,
                          cmd="/quobyte/proteomics-grp/fran/engines/pilot_tools/diann/2.7.0/diann-2.7.0/diann-linux",
                          tools_root="/quobyte/proteomics-grp/fran/engines/pilot_tools",
                          platform_class="hpc")
        self.assertIs(rec["mismatch"], True)
        self.assertIn("PIN_ENGINE=diann PIN_VERSION=2.6.1 bash scripts/acquire_tools.sh hpc "
                      "/quobyte/proteomics-grp/fran/engines/pilot_tools", err)

    def _assert_pasteable(self, err, *commands):
        """Each command the message tells the user to run is a line of its own, and bash can
        parse every such line as it stands. The first version glued the explanation onto the
        end -- `... acquire_tools.sh hpc <root> (this rewrites every engine's entry ...)` -- so
        the line an agent or user copied failed with `syntax error near unexpected token '('`.
        Here "a line of its own" means: some stderr line, stripped, IS the command."""
        lines = [ln.strip() for ln in err.splitlines()]
        for cmd in commands:
            self.assertIn(cmd, lines, err)
            r = subprocess.run(["bash", "-n", "-c", cmd], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, f"{cmd!r}: {r.stderr}")

    def test_the_mismatch_fix_it_command_can_be_pasted_as_it_stands(self):
        rec, err = record("diann", {"diann": "2.7.0"}, DIANN_261,
                          cmd="/quobyte/proteomics-grp/fran/engines/pilot_tools/diann/2.7.0/diann-2.7.0/diann-linux",
                          tools_root="/quobyte/proteomics-grp/fran/engines/pilot_tools",
                          platform_class="hpc")
        self.assertIs(rec["mismatch"], True)
        self._assert_pasteable(err, "PIN_ENGINE=diann PIN_VERSION=2.6.1 bash "
                                    "scripts/acquire_tools.sh hpc "
                                    "/quobyte/proteomics-grp/fran/engines/pilot_tools")
        self.assertIn("rewrites every engine's entry", err)      # the caveat is still said

    def test_the_null_note_fix_it_command_can_be_pasted_as_it_stands(self):
        rec, err = record("diann", {"diann": "latest"}, DIANN_261, cmd="/opt/engine/bin",
                          tools_root="/home/u/my tools", platform_class="linux")
        self.assertIsNone(rec["value"])
        self._assert_pasteable(err, "PIN_ENGINE=diann PIN_VERSION=2.6.1 bash "
                                    "scripts/acquire_tools.sh linux '/home/u/my tools'")

    def test_the_mac_fix_it_commands_can_be_pasted_as_they_stand(self):
        rec, err = record("diann", {"diann": "2.7.0"}, DIANN_261,
                          cmd="docker run --rm -v $PWD:/data proteomics-pipeline/diann:2.7.0 diann-linux",
                          platform_class="mac", tools_root="/Users/u/.proteomics-pipeline/tools")
        self.assertIs(rec["mismatch"], True)
        self._assert_pasteable(err, "bash scripts/build_diann_docker.sh 2.6.1",
                               "DIANN_DOCKER_IMAGE=proteomics-pipeline/diann:2.6.1 bash "
                               "scripts/acquire_tools.sh mac /Users/u/.proteomics-pipeline/tools")

    def test_mac_docker_mismatch_fix_it_rebuilds_the_image(self):
        """On macOS acquire_tools.sh only wraps $DIANN_DOCKER_IMAGE, so re-running it cannot
        change the version; the image has to be built for the pin."""
        rec, err = record("diann", {"diann": "2.7.0"}, DIANN_261,
                          cmd="docker run --rm -v $PWD:/data proteomics-pipeline/diann:2.7.0 diann-linux",
                          platform_class="mac", tools_root="/Users/u/.proteomics-pipeline/tools")
        self.assertIs(rec["mismatch"], True)
        self.assertIn("build_diann_docker.sh 2.6.1", err)

    def test_another_engines_pin_is_never_borrowed(self):
        """The manifest pins ONE engine. tools.json has no `versions` entry for AlphaDIA or
        FragPipe, and `--engine alphadia` under a DIA-NN manifest used to be stamped with
        DIA-NN's version."""
        rec, _ = record("alphadia", {"diann": "2.6.1", "sage": "latest"},
                        {"name": "diann", "version": "2.6.1"},
                        cmd="apptainer exec --bind /quobyte:/quobyte "
                            "/quobyte/proteomics-grp/apptainers/alphadia.sif alphadia")
        self.assertIsNone(rec["value"])
        self.assertIsNone(rec["manifest_pin"])
        self.assertIsNone(rec["mismatch"])

    def test_no_version_anywhere_is_null_and_said(self):
        rec, err = record("fragpipe", {}, {"name": "fragpipe"})
        self.assertIsNone(rec["value"])
        self.assertTrue(rec["source"].startswith("unknown -- "), rec["source"])
        self.assertIn("fragpipe", err)

    # -- what counts as a version at all -----------------------------------------------
    # The first version of this check listed the strings acquire_tools.sh writes today
    # ("latest", "env", "unknown", ""). Everything else went through: `versions.sage:
    # "nightly"` was recorded as the engine version of a deposited search. tools.json is
    # written by a shell script whose every "cannot tell" branch is one edit away from a new
    # word, so the test is the SHAPE of a version, not a list of the known non-versions.
    NOT_VERSIONS = ("latest", "env", "unknown", "nightly", "dev", "n/a", "null", "TBD",
                    "None", "<unknown>", "garbage", "main", "HEAD", "2.x", "v", "2026-09-16")

    def test_a_string_that_is_not_shaped_like_a_version_is_never_the_version(self):
        for bad in self.NOT_VERSIONS:
            rec, err = record("diann", {"diann": bad}, DIANN_261)
            self.assertIsNone(rec["value"], bad)
            self.assertTrue(rec["source"].startswith("unknown -- "), (bad, rec["source"]))
            self.assertIn(repr(bad), rec["source"], bad)     # the string itself is on the record
            self.assertNotIn("2.6.1", rec["source"], bad)    # and never the pin
            self.assertIsNone(rec["mismatch"], bad)
            self.assertEqual(rec["tools_json"], bad, bad)
            self.assertIn("null", err, bad)

    def test_a_version_shape_is_still_read_as_one(self):
        """The refusal above must not swallow the answers: a release number, with or without
        the `v` Sage's tags carry, and however many components it has."""
        for good, want in (("2.6.1", "2.6.1"), ("v0.14.7", "0.14.7"), ("24.0", "24.0"),
                           ("1.2.3.4", "1.2.3.4"), (" 2.6.1 ", "2.6.1")):
            rec, _ = record("diann", {"diann": good}, None)
            self.assertEqual(rec["value"], want, good)
            self.assertTrue(rec["source"].startswith("tools.json"), (good, rec["source"]))

    def test_a_non_version_in_tools_json_still_lets_the_command_name_the_build(self):
        """"nightly" is no more a build than "latest" is: the HIVE path the command runs
        still names 2.6.0, and that is what ran."""
        rec, _ = record("diann", {"diann": "nightly"}, DIANN_261,
                        cmd="/quobyte/proteomics-grp/dia-nn/build_260/diann-2.6.0/diann-linux")
        self.assertEqual(rec["value"], "2.6.0")
        self.assertTrue(rec["source"].startswith("named by the command"), rec["source"])

    def test_a_nameless_manifest_engine_block_pins_nothing(self):
        """`{"engine": {"version": "2.6.1"}}` names no engine, so it is not a pin for THIS
        one. Read as one, it stamped a Sage search with `mismatch: true` against DIA-NN's
        number and warned that "the manifest pins sage 2.6.1"."""
        rec, err = record("sage", {"sage": ""}, {"version": "2.6.1"},
                          cmd="/opt/sage-v0.14.7-x86_64-unknown-linux-gnu/sage")
        self.assertEqual(rec["value"], "0.14.7")
        self.assertIsNone(rec["manifest_pin"])
        self.assertIsNone(rec["mismatch"])
        self.assertNotIn("WARNING", err)
        self.assertNotIn("2.6.1", err)

    def test_a_command_that_names_two_builds_is_refused_not_resolved_to_the_last_one(self):
        """An apptainer line names the image AND the binary inside it. When they disagree,
        nothing here can tell which one ran -- the refusal below is what that case is for,
        and keeping only the last match per string hid it."""
        rec, err = record("diann", {"diann": ""}, DIANN_261,
                          cmd="apptainer exec --bind /quobyte:/quobyte "
                              "/quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif "
                              "/diann-2.6.1/diann-linux")
        self.assertEqual(rec["named_by_command"], ["2.3.0", "2.6.1"])
        self.assertIsNone(rec["value"])
        self.assertIsNone(rec["mismatch"])
        self.assertIn("several", rec["source"])
        self.assertIn("null", err)

    def test_nothing_but_a_version_number_can_ever_become_the_value(self):
        """The property the whole record exists for, over every combination of junk the three
        inputs can hold: whatever `value` is, it is a version number or it is None. It is the
        one field fran_deposit.py forwards to FRAN, where it is stored as the engine version
        of a deposited search and cited in published results."""
        junk = ("latest", "env", "unknown", "nightly", "dev", "n/a", "null", "TBD", "None",
                "<unknown>", "garbage", "", "  ", "v", "2.x", "2026-09-16", "0x2", "2",
                "sage 0.14.6", "2.6.1-rc1", "diann", "../../etc/passwd", "2.6.1; rm -rf /")
        cmds = ("/opt/engine/bin", "/opt/dia-nn/nightly/diann-linux", "docker run diann:latest",
                "/opt/dia-nn/build_260/diann-2.6.0/diann-linux")
        shape = re.compile(r"\d+(\.\d+)+")
        for tj in junk:
            for pin in junk:
                for cmd in cmds:
                    rec, _ = record("diann", {"diann": tj}, {"name": "diann", "version": pin},
                                    cmd=cmd)
                    v = rec["value"]
                    self.assertTrue(v is None or shape.fullmatch(v), (tj, pin, cmd, v))
                    self.assertNotEqual(v, pin, (tj, pin, cmd))   # never the bare request

    def test_the_record_has_exactly_the_provenance_record_shape(self):
        """One object with `value` and `source`, like scan_window -- not a second, flat set of
        top-level keys beside main's records."""
        rec, _ = record("diann", {"diann": "2.7.0"}, DIANN_261)
        self.assertEqual(set(rec), {"value", "source", "tools_json", "named_by_command",
                                    "manifest_pin", "mismatch"})


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
            ev = prov["engine_version"]
            self.assertEqual(ev["value"], "2.7.0")
            self.assertEqual(ev["tools_json"], "2.7.0")
            self.assertEqual(ev["manifest_pin"], "2.6.1")
            self.assertIs(ev["mismatch"], True)
            self.assertTrue(ev["source"].startswith("tools.json"), ev["source"])
            self.assertIn("WARNING", err)

    def test_matching_versions_record_no_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            prov, _ = self._run(d, "2.6.1", "2.6.1")
            self.assertIs(prov["engine_version"]["mismatch"], False)
            self.assertEqual(prov["version"], "2.6.1")

    def test_one_record_beside_the_others_not_a_parallel_set_of_keys(self):
        """search_provenance.json describes each thing with one self-describing object
        (`scan_window`, `sbatch_refused`). The engine version is one more such object; the flat
        `version_source` / `*_engine_version` / `engine_version_mismatch` keys an earlier
        revision of this change spread over the top level must not come back beside it."""
        with tempfile.TemporaryDirectory() as d:
            prov, _ = self._run(d, "latest", "2.6.1",
                                diann_cmd="/quobyte/proteomics-grp/dia-nn/build_260/diann-2.6.0/diann-linux")
            self.assertEqual(prov["version"], prov["engine_version"]["value"])
            self.assertEqual(prov["version"], "2.6.0")
            self.assertIn("scan_window", prov)
            for flat in ("version_source", "tools_engine_version", "command_engine_versions",
                         "manifest_engine_version", "engine_version_mismatch"):
                self.assertNotIn(flat, prov)

    def test_the_chain_route_records_the_engine_version_too(self):
        """A DIA-NN search of >5 files on a SLURM host routes to the 5-step chain, and with
        --sbatch run_search.py exits 3 after writing search_provenance.json. The engine
        version must be in THAT file as well -- most real cohorts take this route."""
        with tempfile.TemporaryDirectory() as d:
            bindir = os.path.join(d, "bin")
            os.makedirs(bindir)
            with open(os.path.join(bindir, "sbatch"), "w") as fh:      # slurm_available()
                fh.write("#!/bin/sh\necho 1\n")
            os.chmod(os.path.join(bindir, "sbatch"), 0o755)
            raws = []
            for i in range(6):
                raws.append(os.path.join(d, f"f{i}.d"))
                os.makedirs(raws[-1])
            cfg = os.path.join(d, "diann.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n--xic 10\n--mobilograms\n")
            fasta = os.path.join(d, "db.fasta")
            with open(fasta, "w") as fh:
                fh.write(">sp|P1|X\nPEPTIDEK\n")
            tools = os.path.join(d, "tools.json")
            with open(tools, "w") as fh:
                json.dump({"diann": "/opt/dia-nn/build_270/diann-2.7.0/diann-linux",
                           "versions": {"diann": "2.7.0"}}, fh)
            bundle = os.path.join(d, "workflow.manifest.json")
            with open(bundle, "w") as fh:
                json.dump({"acquisition": "DIA",
                           "engine": {"name": "diann", "version": "2.6.1"}}, fh)
            out = os.path.join(d, "search_out")
            env = {k: v for k, v in os.environ.items() if k != "SLURM_JOB_ID"}
            env["PATH"] = bindir + os.pathsep + "/usr/bin:/bin"
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                                "--tools", tools, "--bundle", bundle, "--params", cfg,
                                "--fasta", fasta, "--out", out, "--files", *raws,
                                "--sbatch", "job.sh"],
                               cwd=d, capture_output=True, text=True, env=env, timeout=120)
            self.assertEqual(r.returncode, run_search.SBATCH_NOT_WRITTEN, r.stdout + r.stderr)
            with open(os.path.join(out, "search_provenance.json")) as fh:
                prov = json.load(fh)
            self.assertEqual(prov["search_mode"], "parallel_5step")
            self.assertEqual(prov["version"], "2.7.0")
            self.assertIs(prov["engine_version"]["mismatch"], True)
            self.assertEqual(prov["engine_version"]["manifest_pin"], "2.6.1")
            self.assertIn("step 1b", prov["scan_window"]["source"])
            self.assertIn("engine version mismatch", r.stderr)

    def test_a_non_version_in_tools_json_reaches_the_file_as_null(self):
        """End to end: the string acquire_tools.sh wrote is kept in the record, and the
        `version` FRAN is given is null -- not the word."""
        with tempfile.TemporaryDirectory() as d:
            prov, err = self._run(d, "nightly", "2.6.1")
            self.assertIsNone(prov["version"])
            ev = prov["engine_version"]
            self.assertIsNone(ev["value"])
            self.assertEqual(ev["tools_json"], "nightly")
            self.assertEqual(ev["manifest_pin"], "2.6.1")
            self.assertIsNone(ev["mismatch"])
            self.assertIn("not a version", ev["source"])
            self.assertIn("null", err)

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
