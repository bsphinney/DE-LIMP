#!/usr/bin/env python3
"""
Finding ThermoRawFileParser, and the .NET a framework-dependent build of it runs on.

gabrig 2026-09-23, HIVE, 15 Fusion Lumos .raw, a Core member:
  - detect_acquisition.py said "ThermoRawFileParser not found" for all 15, although the Core
    keeps TRFP 2.0.0.0 at /quobyte/proteomics-grp/tools/ThermoRawFileParser/ -- nothing
    looked there, and setup.json still said ready_for.dia. The search would have run the
    380-980 FALLBACK on a method that acquired 357-1105;
  - that copy is framework-dependent: its runtimeconfig lists Microsoft.NETCore.App 8.0.0 AND
    Microsoft.AspNetCore.App 8.0.0. With no .NET it exits 131, and with ensure_dotnet8.sh's
    NETCore-only install as DOTNET_ROOT it exits 150 "You must install or update .NET" --
    once per file, after two srun jobs;
  - an Orbitrap .raw carries no resolution, so estimate_params.py classes it orbitrap_generic
    and the parallel chain declines -- and nothing told the orchestrator to ask.

No real parser, no real .NET and no network: the parser is a shell shim that behaves like the
.NET host (131 with no root, 150 without AspNetCore) and then replays the captured TRFP
output (fixtures/trfp/fake_trfp.py); a .NET root is a directory tree
(host/fxr/<v>, shared/<framework>/<v>) with a `dotnet` that lists it; `curl` is a stub that
hands ensure_dotnet8.sh a stub dotnet-install.sh. POSIX only.
"""
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import detect_acquisition as da           # noqa: E402
import estimate_params                     # noqa: E402
# Module imports, not `from ... import`: a TestCase class bound in this namespace would be
# collected and run a second time here.
import test_thermo_raw_detection as trfp   # noqa: E402
import test_setup_diann_readiness as setup_t  # noqa: E402

NETCORE, ASPNET = "Microsoft.NETCore.App", "Microsoft.AspNetCore.App"
TRFP2_RUNTIMECONFIG = {"runtimeOptions": {"tfm": "net8.0", "frameworks": [
    {"name": NETCORE, "version": "8.0.0"}, {"name": ASPNET, "version": "8.0.0"}]}}


def _exe(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


# `dotnet --list-runtimes` read off the directory tree it sits in, so the list and the
# directories can never disagree. Pure bash (bash 3.2), no coreutils.
FAKE_DOTNET = r'''#!/bin/bash
root="$(cd "${0%/*}" && pwd)"
[ "$1" = --list-runtimes ] || exit 0
for d in "$root"/shared/*/*/; do
  [ -d "$d" ] || continue
  d="${d%/}"; v="${d##*/}"; n="${d%/*}"; n="${n##*/}"
  echo "$n $v [$root/shared/$n]"
done
'''


def fake_dotnet_root(path, frameworks=(NETCORE, ASPNET), version="8.0.28"):
    """A .NET install as the host sees it: host/fxr/<v> and shared/<framework>/<v>."""
    os.makedirs(os.path.join(path, "host", "fxr", version), exist_ok=True)
    for fw in frameworks:
        os.makedirs(os.path.join(path, "shared", fw, version), exist_ok=True)
    _exe(os.path.join(path, "dotnet"), FAKE_DOTNET)
    return path


# The TRFP apphost as the .NET host runs it: no usable DOTNET_ROOT -> 131, a root without
# every framework -> 150 (the messages are the host's own), else the replaying stand-in.
# Every start is logged before the host checks, so a test can tell "never started" apart
# from "started and refused".
def fake_apphost(path, requires=(NETCORE, ASPNET)):
    checks = "\n".join(
        f'has "{fw}" || {{ echo "You must install or update .NET to run this application." >&2; '
        f'exit 150; }}' for fw in requires)
    return _exe(path, f'''#!/bin/sh
[ -n "$SHIM_LOG" ] && echo "$*" >> "$SHIM_LOG"
r="$DOTNET_ROOT"
if [ -z "$r" ] || [ ! -d "$r/host/fxr" ]; then
  echo "You must install .NET to run this application." >&2; exit 131
fi
has() {{ for d in "$r"/shared/"$1"/8.*; do [ -d "$d" ] && return 0; done; return 1; }}
{checks}
exec "{sys.executable}" "{trfp.FAKE}" "$@"
''')


def _lines(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return fh.read().splitlines()


def write_runtimeconfig(apphost, config):
    with open(apphost + ".runtimeconfig.json", "w") as fh:
        json.dump(config, fh)


class _DotnetCase(trfp._FakeParserCase):
    """_FakeParserCase with no parser on PATH, no .NET anywhere the machine running the tests
    could supply one, and a framework-dependent TRFP 2.x apphost ready to be put somewhere."""

    def setUp(self):
        super().setUp()
        self.hide_parser()
        os.environ["PATH"] = os.pathsep.join(
            d for d in os.environ["PATH"].split(os.pathsep)
            if not os.path.exists(os.path.join(d, "dotnet")))
        for k in ("DOTNET_ROOT", "DOTNET_CORE_SDK_ROOT", "DOTNET_ROOT_X64", "DOTNET_ROOT_ARM64"):
            os.environ.pop(k, None)
        os.environ["PROTEOMICS_DOTNET_DIR"] = os.path.join(self.tmp, "pp_dotnet8")
        os.environ["PROTEOMICS_DOTNET_SYSTEM_ROOTS"] = ""
        self.shim_log = os.path.join(self.tmp, "shim_starts.log")
        os.environ["SHIM_LOG"] = self.shim_log
        self.core = os.path.join(self.tmp, "core", "ThermoRawFileParser")
        fake_apphost(self.core)
        write_runtimeconfig(self.core, TRFP2_RUNTIMECONFIG)

    def starts(self):
        if not os.path.exists(self.shim_log):
            return []
        with open(self.shim_log) as fh:
            return [ln.rstrip("\n") for ln in fh]


# --------------------------------------------------------------------------------------------
class ParserIsFoundWhereItIs(_DotnetCase):
    def test_the_core_shared_copy_is_found_last(self):
        os.environ["THERMORAWFILEPARSER_SHARED"] = self.core
        self.assertEqual(da.locate_trfp()[0], [self.core])
        self.assertIn("THERMORAWFILEPARSER_SHARED", da.locate_trfp()[1])

    def test_the_default_shared_copy_is_the_core_one(self):
        """The one path gabrig's session needed, and nothing looked at."""
        self.assertEqual(da.TRFP_SHARED_DEFAULT,
                         "/quobyte/proteomics-grp/tools/ThermoRawFileParser/ThermoRawFileParser")
        os.environ.pop("THERMORAWFILEPARSER_SHARED")
        self.assertEqual(da._shared_trfp_paths(), [da.TRFP_SHARED_DEFAULT])

    def test_a_shared_list_is_searched_in_order_and_skips_what_is_not_there(self):
        os.environ["THERMORAWFILEPARSER_SHARED"] = os.pathsep.join(
            [os.path.join(self.tmp, "nope", "ThermoRawFileParser"), self.core])
        self.assertEqual(da.find_trfp(), [self.core])

    def test_order_is_env_var_then_path_then_pipeline_env_then_shared(self):
        os.environ["THERMORAWFILEPARSER_SHARED"] = self.core
        envbin = os.path.join(self.tmp, "env", "bin")
        in_env = _exe(os.path.join(envbin, "ThermoRawFileParser"), "#!/bin/sh\n")
        with open(os.path.join(self.tmp, "setup.json"), "w") as fh:   # PROTEOMICS_PIPELINE_HOME
            json.dump({"env_prefix": os.path.dirname(envbin)}, fh)
        self.assertEqual(da.locate_trfp(), ([in_env], "the pipeline's conda env"))

        on_path = _exe(os.path.join(self.tmp, "pathbin", "ThermoRawFileParser"), "#!/bin/sh\n")
        os.environ["PATH"] = os.path.dirname(on_path) + os.pathsep + os.environ["PATH"]
        self.assertEqual(da.locate_trfp(), ([on_path], "PATH"))

        os.environ["THERMORAWFILEPARSER"] = "dotnet /opt/trfp/ThermoRawFileParser.dll"
        self.assertEqual(da.locate_trfp()[0], ["dotnet", "/opt/trfp/ThermoRawFileParser.dll"])

    def test_not_found_names_the_shared_copy_it_looked_for_and_the_public_sources(self):
        missing = os.path.join(self.tmp, "nope", "ThermoRawFileParser")
        os.environ["THERMORAWFILEPARSER_SHARED"] = missing
        r = da.classify(self.raw(trfp.LUMOS))
        self.assertIn("not found", r["reason"])
        self.assertIn(missing, r["reason"])
        self.assertIn("github.com/compomics/ThermoRawFileParser", r["reason"])
        self.assertIn("bioconda thermorawfileparser", r["reason"])
        self.assertIn("setup.sh", r["reason"])


class DotnetRootIsChosenByFramework(_DotnetCase):
    def setUp(self):
        super().setUp()
        os.environ["THERMORAWFILEPARSER_SHARED"] = self.core
        self.netcore_only = fake_dotnet_root(os.path.join(self.tmp, "netcore_only"), (NETCORE,))
        self.both = fake_dotnet_root(os.path.join(self.tmp, "both"))

    def test_the_needs_are_read_from_the_runtimeconfig(self):
        needs, muxer, rc = da.dotnet_frameworks_needed([self.core])
        self.assertEqual([n[0] for n in needs], [NETCORE, ASPNET])
        self.assertEqual(needs[0][1], (8, 0, 0))
        self.assertFalse(muxer)
        self.assertEqual(rc, os.path.realpath(self.core) + ".runtimeconfig.json")

    def test_a_netcore_only_root_is_rejected_and_one_with_both_is_used(self):
        """ensure_dotnet8.sh's old install as $DOTNET_ROOT is exactly the exit-150 case."""
        os.environ["DOTNET_ROOT"] = self.netcore_only
        os.environ["PROTEOMICS_DOTNET_DIR"] = self.both
        launch = da.trfp_launch([self.core])
        self.assertIsNone(launch["problem"])
        self.assertEqual(launch["dotnet_root"], self.both)
        self.assertEqual(da.dotnet_root_lacks(self.netcore_only, da.dotnet_frameworks_needed(
            [self.core])[0]), [f"no {ASPNET} 8.0.0+"])

    def test_the_hive_module_root_is_picked_up(self):
        """`module load dotnet-core-sdk/8.0.4` sets DOTNET_CORE_SDK_ROOT, not DOTNET_ROOT; its
        8.0.4 has both frameworks, and 8.0.4 is fine for TRFP (only DIA-NN wants >= 8.0.17)."""
        module = fake_dotnet_root(os.path.join(self.tmp, "spack", "dotnet-core-sdk-8.0.4"),
                                  version="8.0.4")
        os.environ["DOTNET_ROOT"] = self.netcore_only
        os.environ["DOTNET_CORE_SDK_ROOT"] = module
        self.assertEqual(da.trfp_launch([self.core])["dotnet_root"], module)

    def test_dotnet_on_path_is_followed_to_its_real_root(self):
        link = os.path.join(self.tmp, "usrbin", "dotnet")
        os.makedirs(os.path.dirname(link))
        os.symlink(os.path.join(self.both, "dotnet"), link)
        os.environ["PATH"] = os.path.dirname(link) + os.pathsep + os.environ["PATH"]
        self.assertEqual(da.trfp_launch([self.core])["dotnet_root"], os.path.realpath(self.both))

    def test_a_newer_major_does_not_count(self):
        """rollForward defaults to Minor: a .NET 9 install does not start a net8.0 app."""
        nine = fake_dotnet_root(os.path.join(self.tmp, "nine"), version="9.0.10")
        os.environ["DOTNET_ROOT"] = nine
        launch = da.trfp_launch([self.core])
        self.assertIsNotNone(launch["problem"])
        self.assertIn("(has 9.0.10)", launch["problem"])

    def test_the_chosen_root_reaches_the_parser_and_the_reader_record(self):
        os.environ["DOTNET_ROOT"] = self.netcore_only          # would be exit 150 as is
        os.environ["PROTEOMICS_DOTNET_DIR"] = self.both
        r = da.classify(self.raw(trfp.LUMOS))
        self.assertEqual(r["acquisition"], "DIA", r["reason"])
        self.assertEqual(r["instrument"], "Orbitrap Fusion Lumos")
        self.assertEqual(r["precursor_mz_range"], [350.05, 1200.95])
        self.assertTrue(r["reader"].startswith("ThermoRawFileParser 2.0.0.0 "), r["reader"])
        self.assertIn(f"DOTNET_ROOT={self.both}", r["reader"])
        self.assertEqual(r["warnings"], [])
        self.assertEqual(os.environ["DOTNET_ROOT"], self.netcore_only,
                         "the parser's DOTNET_ROOT leaked into this process")

    def test_a_self_contained_build_needs_no_dotnet(self):
        write_runtimeconfig(self.core, {"runtimeOptions": {"tfm": "net8.0", "includedFrameworks": [
            {"name": NETCORE, "version": "8.0.21"}]}})
        launch = da.trfp_launch([self.core])
        self.assertEqual((launch["needs"], launch["dotnet_root"], launch["problem"]),
                         ([], None, None))

    def test_dotnet_dll_form_runs_under_a_dotnet_that_has_both(self):
        """The muxer resolves frameworks from its OWN directory and ignores DOTNET_ROOT, so the
        dotnet that is named is replaced by the one whose root qualifies."""
        dll = os.path.join(self.tmp, "net8zip", "ThermoRawFileParser.dll")
        os.makedirs(os.path.dirname(dll))
        open(dll, "w").close()
        write_runtimeconfig(dll[:-4], TRFP2_RUNTIMECONFIG)
        cmd = [os.path.join(self.netcore_only, "dotnet"), dll]
        needs, muxer, _ = da.dotnet_frameworks_needed(cmd)
        self.assertTrue(muxer)
        self.assertEqual(len(needs), 2)
        os.environ["PROTEOMICS_DOTNET_DIR"] = self.both
        launch = da.trfp_launch(cmd)
        self.assertEqual(launch["cmd"], [os.path.join(self.both, "dotnet"), dll])
        self.assertEqual(launch["dotnet_root"], self.both)


class ACannotStartParserFailsOnceUpFront(_DotnetCase):
    def setUp(self):
        super().setUp()
        os.environ["THERMORAWFILEPARSER_SHARED"] = self.core
        os.environ["DOTNET_ROOT"] = fake_dotnet_root(os.path.join(self.tmp, "nc"), (NETCORE,))

    def _run(self, *args):
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                               *args], capture_output=True, text=True, env=os.environ.copy())

    def test_one_message_with_the_fix_and_no_parser_started(self):
        raws = [self.raw(trfp.LUMOS, subdir=f"b{i}") for i in range(3)]
        res = self._run(*raws)
        self.assertEqual(res.returncode, 0, res.stderr)
        errors = [ln for ln in res.stderr.splitlines() if "ERROR:" in ln]
        self.assertEqual(len(errors), 1, res.stderr)
        msg = errors[0]
        for words in (ASPNET, "ensure_dotnet8.sh", "module load dotnet-core-sdk/8.0.4",
                      os.environ["DOTNET_ROOT"], "FALLBACK", "exits 150"):
            self.assertIn(words, msg)
        self.assertEqual(self.starts(), [], "a parser that cannot start was started anyway")
        payload = json.loads(res.stdout)
        self.assertTrue(payload["needs_confirmation"])
        reasons = {f["reason"] for f in payload["files"]}
        self.assertEqual(len(reasons), 1, "every file must say the same thing")
        self.assertIn(reasons.pop(), msg)
        for f in payload["files"]:
            self.assertEqual(f["acquisition"], "unknown")
            self.assertEqual(f["warnings"], [f["reason"]])
            self.assertIn("cannot start", f["reader"])
        # the per-file stderr lines point at the one message instead of repeating it
        warn = [ln for ln in res.stderr.splitlines() if "WARNING:" in ln]
        self.assertEqual(len(warn), 3, res.stderr)
        self.assertTrue(all(da.CANNOT_START_SEE_ABOVE in ln and ASPNET not in ln
                            for ln in warn), warn)

    def test_on_a_login_node_it_is_not_sent_to_srun_first(self):
        """gabrig's two srun jobs: the refusal sent a cohort to a compute node only for every
        file to fail there the way a directory check on the login node already knew."""
        sbin = os.path.join(self.tmp, "sbin")
        _exe(os.path.join(sbin, "sbatch"), "#!/bin/sh\nexit 0\n")
        os.environ["PATH"] = sbin + os.pathsep + os.environ["PATH"]
        os.environ.pop("SLURM_JOB_ID", None)
        self.assertTrue(da.on_cluster_login_node())
        res = self._run(*[self.raw(trfp.LUMOS, subdir=f"b{i}")
                          for i in range(da.LOGIN_NODE_MAX_RAW + 1)])
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertNotIn("REFUSING", res.stderr)
        self.assertIn(ASPNET, res.stderr)
        self.assertEqual(self.starts(), [])

    def test_a_host_refusal_the_directories_did_not_predict_is_caught_by_one_version_call(self):
        """A runtimeconfig that under-states the needs (or a host quirk): the directory check
        passes, the host still exits 150. One `--version`, then the same single message."""
        write_runtimeconfig(self.core, {"runtimeOptions": {"tfm": "net8.0", "frameworks": [
            {"name": NETCORE, "version": "8.0.0"}]}})
        res = self._run(*[self.raw(trfp.LUMOS, subdir=f"b{i}") for i in range(2)])
        self.assertEqual(res.returncode, 0, res.stderr)
        errors = [ln for ln in res.stderr.splitlines() if "ERROR:" in ln]
        self.assertEqual(len(errors), 1, res.stderr)
        self.assertIn("exit 150", errors[0])
        self.assertIn("You must install or update .NET", errors[0])
        self.assertIn("ensure_dotnet8.sh", errors[0])
        self.assertEqual(self.starts(), ["--version"], "only the one start check may run")

    def test_it_still_reads_once_the_fix_is_in(self):
        os.environ["PROTEOMICS_DOTNET_DIR"] = fake_dotnet_root(os.path.join(self.tmp, "fixed"))
        res = self._run(self.raw(trfp.LUMOS))
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertNotIn("ERROR:", res.stderr)
        f = json.loads(res.stdout)["files"][0]
        self.assertEqual((f["acquisition"], f["instrument"]), ("DIA", "Orbitrap Fusion Lumos"))


class ReaderCheckForSetup(_DotnetCase):
    def _check(self):
        res = subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                              "--check-reader"], capture_output=True, text=True,
                             env=os.environ.copy())
        return res.returncode, json.loads(res.stdout)

    def test_not_found_is_not_ready(self):
        rc, st = self._check()
        self.assertEqual((rc, st["ready"], st["command"]), (1, False, None))
        self.assertIn("not found", st["note"])

    def test_missing_aspnetcore_is_not_ready_and_says_the_fix(self):
        os.environ["THERMORAWFILEPARSER_SHARED"] = self.core
        rc, st = self._check()
        self.assertEqual((rc, st["ready"]), (1, False))
        self.assertEqual(st["command"], self.core)
        self.assertEqual(st["dotnet_needs"], [f"{NETCORE} 8.0.0+", f"{ASPNET} 8.0.0+"])
        self.assertIn("ensure_dotnet8.sh", st["note"])
        self.assertEqual(self.starts(), [])

    def test_ready_names_the_root_and_the_version(self):
        os.environ["THERMORAWFILEPARSER_SHARED"] = self.core
        root = fake_dotnet_root(os.environ["PROTEOMICS_DOTNET_DIR"])
        rc, st = self._check()
        self.assertEqual((rc, st["ready"], st["version"], st["dotnet_root"]),
                         (0, True, "2.0.0.0", root))
        self.assertEqual(self.starts(), ["--version"], "readiness must not read any .raw")


class OrbitrapResolutionMustBeAsked(trfp._FakeParserCase):
    def _run(self, *paths):
        res = subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                              *paths], capture_output=True, text=True, env=os.environ.copy())
        self.assertEqual(res.returncode, 0, res.stderr)
        return json.loads(res.stdout)

    def test_orbitrap_raws_are_listed_with_the_question_to_ask(self):
        lumos, exploris = self.raw(trfp.LUMOS), self.raw(trfp.EXPLORIS)
        out = self._run(lumos, exploris)
        hint = out["orbitrap_resolution_unknown"]
        self.assertEqual(sorted(hint["files"]), sorted([lumos, exploris]))
        self.assertEqual(hint["instruments"], ["Orbitrap Exploris 480", "Orbitrap Fusion Lumos"])
        for words in ("ASK the user", "--ms1-resolution", "--ms2-resolution",
                      "orbitrap_generic", "parallel chain declines"):
            self.assertIn(words, hint["ask"])

    def test_it_is_not_a_warning_and_does_not_ask_for_confirmation(self):
        """A clean read stays clean: acquisition and range stand whether or not the resolution
        came with them, and a warning would put a confirmation on every such Orbitrap cohort."""
        out = self._run(self.raw(trfp.EXPLORIS))
        self.assertIsNotNone(out["orbitrap_resolution_unknown"])
        self.assertFalse(out["needs_confirmation"])
        self.assertEqual(out["files"][0]["warnings"], [])

    def test_nothing_to_ask_without_an_orbitrap_raw(self):
        mzml = os.path.join(self.tmp, "run.mzML")
        open(mzml, "w").close()
        self.assertIsNone(self._run(mzml)["orbitrap_resolution_unknown"])

    def test_the_orbitrap_rule_is_estimate_params_own(self):
        """One keyword list (DE-LIMP rule 3): the class is the one estimate_params.py acts on."""
        self.assertIs(da.classify_instrument, estimate_params.classify_instrument)
        self.assertEqual(estimate_params.classify_instrument("Orbitrap Fusion Lumos")[0],
                         "orbitrap_generic")
        self.assertIsNone(da.resolution_summary([
            {"vendor": "Thermo", "instrument": "Orbitrap Astral", "file": "a.raw",
             "ms1_resolution": None, "ms2_resolution": None}])["orbitrap_resolution_unknown"],
            "Astral has a documented tolerance: nothing to ask")


# --------------------------------------------------------------------------------------------
ENSURE = os.path.join(SCRIPTS, "ensure_dotnet8.sh")
# `curl -sSL <url> -o <file>` -> <file> is the stub installer below; logs "curl"
FAKE_CURL = '''#!/bin/bash
out=""; while [ $# -gt 0 ]; do [ "$1" = -o ] && out="$2"; shift; done
echo curl >> "$FAKE_LOG"
cp "$FAKE_INSTALLER" "$out"
'''
# dotnet-install.sh --channel 8.0 --runtime <dotnet|aspnetcore> --install-dir <dir>. Like the
# real aspnetcore runtime archive it brings the matching Microsoft.NETCore.App along.
FAKE_INSTALLER = '''#!/bin/bash
rt=""; dir=""
while [ $# -gt 0 ]; do
  case "$1" in --runtime) rt="$2"; shift ;; --install-dir) dir="$2"; shift ;; esac; shift
done
echo "install $rt" >> "$FAKE_LOG"
case ",$FAKE_INSTALL_FAIL," in *",$rt,"*) echo "download failed" >&2; exit 1 ;; esac
v="${FAKE_DOTNET_VERSION:-8.0.28}"
mkdir -p "$dir/host/fxr/$v" "$dir/shared/Microsoft.NETCore.App/$v"
[ "$rt" = aspnetcore ] && mkdir -p "$dir/shared/Microsoft.AspNetCore.App/$v"
cp "$FAKE_DOTNET" "$dir/dotnet"; chmod +x "$dir/dotnet"
'''


class EnsureDotnet8(unittest.TestCase):
    """ensure_dotnet8.sh with a stub curl and a stub installer: no network, ever."""

    TOOLS = ("awk", "bash", "chmod", "cp", "dirname", "mkdir", "mktemp", "readlink", "rm")

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.bin = os.path.join(self.d, "bin")
        os.makedirs(self.bin)
        for tool in self.TOOLS:
            real = shutil.which(tool)
            if not real:
                self.skipTest(f"{tool} is not installed")
            os.symlink(real, os.path.join(self.bin, tool))
        _exe(os.path.join(self.bin, "curl"), FAKE_CURL)
        self.log = os.path.join(self.d, "log")
        self.dest = os.path.join(self.d, "dotnet8")
        self.env = {"PATH": self.bin, "HOME": self.d, "PROTEOMICS_DOTNET_DIR": self.dest,
                    "FAKE_LOG": self.log,
                    "FAKE_INSTALLER": _exe(os.path.join(self.d, "installer.sh"), FAKE_INSTALLER),
                    "FAKE_DOTNET": _exe(os.path.join(self.d, "fake_dotnet"), FAKE_DOTNET)}

    def tearDown(self):
        self._tmp.cleanup()

    def run_it(self, **env):
        p = subprocess.run(["bash", ENSURE], capture_output=True, text=True,
                           env=dict(self.env, **env), timeout=60)
        return p, _lines(self.log)

    def frameworks(self, root):
        return sorted(os.listdir(os.path.join(root, "shared")))

    def test_a_root_with_both_is_reused_without_a_download(self):
        fake_dotnet_root(self.dest)
        p, log = self.run_it()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.splitlines()[-1], self.dest)
        self.assertEqual(log, [])

    def test_a_netcore_only_install_gets_aspnetcore_added_in_place(self):
        """The state of every ~/.proteomics-pipeline/dotnet8 this script made before."""
        fake_dotnet_root(self.dest, (NETCORE,))
        p, log = self.run_it()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.splitlines()[-1], self.dest)
        self.assertEqual(log, ["curl", "install aspnetcore"], "reinstalled what was there")
        self.assertEqual(self.frameworks(self.dest), [ASPNET, NETCORE])

    def test_an_empty_dest_gets_both(self):
        p, log = self.run_it()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.splitlines()[-1], self.dest)
        self.assertEqual(log, ["curl", "install dotnet", "curl", "install aspnetcore"])
        self.assertEqual(self.frameworks(self.dest), [ASPNET, NETCORE])

    def test_an_old_netcore_is_upgraded_not_accepted(self):
        """DIA-NN's hard >= 8.0.17: a DEST at 8.0.4 must get a new NETCore."""
        fake_dotnet_root(self.dest, version="8.0.4")
        p, log = self.run_it()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("install dotnet", log)

    def test_a_module_dotnet_at_8_0_4_is_not_enough_for_diann(self):
        module = fake_dotnet_root(os.path.join(self.d, "module"), version="8.0.4")
        os.symlink(os.path.join(module, "dotnet"), os.path.join(self.bin, "dotnet"))
        p, log = self.run_it()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.splitlines()[-1], self.dest)
        self.assertIn("install dotnet", log)

    def test_a_system_dotnet_with_both_is_used(self):
        system = fake_dotnet_root(os.path.join(self.d, "system"))
        os.symlink(os.path.join(system, "dotnet"), os.path.join(self.bin, "dotnet"))
        p, log = self.run_it()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.splitlines()[-1], os.path.realpath(system))
        self.assertEqual(log, [])

    def test_no_aspnetcore_download_keeps_diann_working_and_says_so(self):
        """DIA-NN's callers read the last line; a failed AspNetCore download must not take
        away a .raw search that worked before AspNetCore was asked for."""
        fake_dotnet_root(self.dest, (NETCORE,))
        p, log = self.run_it(FAKE_INSTALL_FAIL="aspnetcore")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.splitlines()[-1], self.dest)
        self.assertIn("WARNING", p.stderr)
        self.assertIn("AspNetCore", p.stderr)
        self.assertIn("exit 150", p.stderr)

    def test_nothing_usable_is_an_error(self):
        p, _log = self.run_it(FAKE_INSTALL_FAIL="dotnet,aspnetcore")
        self.assertEqual(p.returncode, 1)
        self.assertEqual(p.stdout, "")
        self.assertIn("ERROR", p.stderr)


# --------------------------------------------------------------------------------------------
class SetupReportsThermoRawReadiness(setup_t.SetupCheckHarness):
    """setup.json `ready_for.thermo_raw` + `thermo_raw_reader`, from `setup.sh --check` on a
    mocked PATH (coreutils + python3 only), PP_HOME and parser locations."""

    def setUp(self):
        super().setUp()
        os.symlink(sys.executable, os.path.join(self.sys, "python3"))
        self.core = os.path.join(self.d, "core", "ThermoRawFileParser")
        fake_apphost(self.core)
        write_runtimeconfig(self.core, TRFP2_RUNTIMECONFIG)
        self.base = {"THERMORAWFILEPARSER_SHARED": self.core,
                     "PROTEOMICS_DOTNET_SYSTEM_ROOTS": ""}

    def test_the_core_copy_without_aspnetcore_is_not_ready_and_the_note_says_the_fix(self):
        s = self.check(extra_env=self.base)
        self.assertIs(s["ready_for"]["thermo_raw"], False)
        r = s["thermo_raw_reader"]
        self.assertIs(r["ready"], False)
        self.assertEqual(r["command"], self.core)
        self.assertIn("THERMORAWFILEPARSER_SHARED", r["source"])
        self.assertIn(ASPNET, r["note"])
        self.assertIn("ensure_dotnet8.sh", r["note"])
        self.assertTrue(any("thermo_raw_reader.note" in n for n in s["notes"]), s["notes"])
        # the keys the orchestrator already gates on are still there
        self.assertEqual(set(s["ready_for"]), {"de", "dia", "dda", "thermo_raw"})

    def test_the_core_copy_with_ensure_dotnet8s_install_is_ready(self):
        root = fake_dotnet_root(os.path.join(self.d, ".proteomics-pipeline", "dotnet8"))
        s = self.check(extra_env=self.base)
        r = s["thermo_raw_reader"]
        self.assertIs(s["ready_for"]["thermo_raw"], True, r)
        self.assertEqual((r["version"], r["dotnet_root"]), ("2.0.0.0", root))

    def test_a_parser_in_the_env_is_found_without_activate(self):
        env_bin = os.path.join(self.home, "micromamba", "envs", "proteomics-pipeline", "bin")
        _exe(os.path.join(env_bin, "ThermoRawFileParser"),
             f'#!/bin/sh\nexec "{sys.executable}" "{trfp.FAKE}" "$@"\n')
        s = self.check(extra_env=self.base)
        r = s["thermo_raw_reader"]
        self.assertIs(s["ready_for"]["thermo_raw"], True, r)
        self.assertEqual(r["command"], os.path.join(env_bin, "ThermoRawFileParser"))
        self.assertIsNone(r["dotnet_root"])

    def test_without_python_it_is_not_ready_and_says_why(self):
        os.unlink(os.path.join(self.sys, "python3"))
        s = self.check(extra_env=self.base)
        self.assertIs(s["ready_for"]["thermo_raw"], False)
        self.assertIn("detect_acquisition.py", s["thermo_raw_reader"]["note"])


# argv logged; `create` makes an env with python + Rscript, `install` puts the parser stand-in
# and/or a pythonnet conda-meta record in it, as asked. Both honour only -p: a named env would
# land wherever a user's .condarc says.
FAKE_MICROMAMBA = r'''#!/bin/bash
echo "$*" >> "$FAKE_MM_LOG"
sub="$1"; shift
prefix=""; pkgs=""
while [ $# -gt 0 ]; do
  case "$1" in -p) prefix="$2"; shift ;; -r|-c) shift ;; -*) ;; *) pkgs="$pkgs $1" ;; esac; shift
done
[ -n "$prefix" ] || { echo "no -p" >&2; exit 2; }
case "$sub" in
  create)  mkdir -p "$prefix/bin" "$prefix/conda-meta"; ln -s "$FAKE_PY" "$prefix/bin/python"
           printf '#!/bin/sh\nexit 0\n' > "$prefix/bin/Rscript"; chmod +x "$prefix/bin/Rscript" ;;
  install) [ -n "${FAKE_MM_INSTALL_FAIL:-}" ] && { echo "nothing provides thermorawfileparser" >&2; exit 1; }
           for p in $pkgs; do case "$p" in
             thermorawfileparser)
               printf '#!/bin/sh\nexec "%s" "%s" "$@"\n' "$FAKE_PY" "$FAKE_TRFP" > "$prefix/bin/ThermoRawFileParser"
               chmod +x "$prefix/bin/ThermoRawFileParser" ;;
             pythonnet) echo '{}' > "$prefix/conda-meta/pythonnet-3.1.0-pyhd8ed1ab_0.json" ;;
             pandas) echo '{}' > "$prefix/conda-meta/pandas-2.3.2-py311h0000000_0.json" ;;
           esac; done ;;
esac
'''


# ensure_dotnet8.sh's contract without its download: messages on stderr, the root on the LAST
# stdout line; FAKE_ENSURE_FAIL makes it fail as it does with no internet.
FAKE_ENSURE = r'''#!/bin/bash
echo call >> "$FAKE_ENSURE_LOG"
if [ -n "${FAKE_ENSURE_FAIL:-}" ]; then
  echo "curl: (6) Could not resolve host: dot.net" >&2
  echo "ERROR: .NET 8 install did not yield a >= 8.0.17 runtime" >&2; exit 1
fi
echo "reusing .NET 8 at $FAKE_ENSURE_ROOT" >&2
echo "$FAKE_ENSURE_ROOT"
'''


class SetupInstallsTheParserSeparately(setup_t.SetupCheckHarness):
    """setup.sh (not --check) with a stub micromamba: the parser and pythonnet are their own
    install, after the env exists, into the env's prefix; only what is missing; its failure
    costs nothing else."""

    def setUp(self):
        super().setUp()
        for tool in ("chmod", "ln"):
            os.symlink(shutil.which(tool), os.path.join(self.sys, tool))
        _exe(os.path.join(self.sys, "micromamba"), FAKE_MICROMAMBA)
        self.ensure = _exe(os.path.join(self.d, "ensure_dotnet8_stub.sh"), FAKE_ENSURE)
        self.ensure_log = os.path.join(self.d, "ensure.log")
        self.dotnet_root = fake_dotnet_root(os.path.join(self.d, ".proteomics-pipeline", "dotnet8"))
        self.mm_log = os.path.join(self.d, "mm.log")
        self.prefix = os.path.join(self.home, "micromamba", "envs", "proteomics-pipeline")

    def setup(self, *args, **extra):
        env = {"PATH": self.sys, "HOME": self.d, "PP_HOME": self.home, "QUOBYTE_DIR": "",
               "THERMORAWFILEPARSER_SHARED": "", "PROTEOMICS_DOTNET_SYSTEM_ROOTS": "",
               "FAKE_MM_LOG": self.mm_log, "FAKE_PY": sys.executable, "FAKE_TRFP": trfp.FAKE,
               "PROTEOMICS_ENSURE_DOTNET8": self.ensure, "FAKE_ENSURE_LOG": self.ensure_log,
               "FAKE_ENSURE_ROOT": self.dotnet_root}
        env.update(extra)
        r = subprocess.run([self.bash, setup_t.SETUP, *args], capture_output=True, text=True,
                           env=env, timeout=300)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(self.home, "setup.json")) as fh:
            s = json.load(fh)
        return s, _lines(self.mm_log)

    def test_installed_on_its_own_into_the_env_prefix_and_then_ready(self):
        s, calls = self.setup()
        creates = [c for c in calls if c.startswith("create")]
        installs = [c for c in calls if c.startswith("install")]
        self.assertEqual(len(creates), 1, calls)
        self.assertNotIn("thermorawfileparser", creates[0], "part of the main solve")
        self.assertEqual(len(installs), 1, calls)
        for c in creates + installs:
            self.assertIn(f"-p {self.prefix}", c)
            self.assertNotIn(" -n ", f" {c} ")
        self.assertIn("-c conda-forge -c bioconda thermorawfileparser pythonnet pandas",
                      installs[0])
        self.assertTrue(os.path.exists(os.path.join(
            self.prefix, "conda-meta", "pythonnet-3.1.0-pyhd8ed1ab_0.json")))
        self.assertIs(s["ready_for"]["thermo_raw"], True, s["thermo_raw_reader"])
        self.assertEqual(s["thermo_raw_reader"]["command"],
                         os.path.join(self.prefix, "bin", "ThermoRawFileParser"))

    def test_a_second_run_installs_nothing(self):
        self.setup()
        _s, calls = self.setup()
        self.assertEqual(len([c for c in calls if c.startswith("install")]), 1, calls)

    def test_only_what_is_missing_is_installed(self):
        """An env from before pythonnet joined the step: the parser is there, pythonnet is not."""
        self.setup()
        os.remove(os.path.join(self.prefix, "conda-meta", "pythonnet-3.1.0-pyhd8ed1ab_0.json"))
        _s, calls = self.setup()
        installs = [c for c in calls if c.startswith("install")]
        self.assertEqual(len(installs), 2, calls)
        self.assertTrue(installs[1].endswith("-c conda-forge -c bioconda pythonnet"), installs[1])

    def test_the_resolution_reader_is_reported(self):
        s, _calls = self.setup(THERMO_RESOLUTION_PYTHON="")
        rr = s["thermo_raw_reader"]["resolution_reader"]
        self.assertIs(rr["ready"], False)
        self.assertIn("pythonnet", rr["note"])
        self.assertIs(s["ready_for"]["thermo_raw"], True, "the parser alone decides that")

    def test_an_existing_env_gets_pandas_on_a_rerun(self):
        """gabrig #15: the HIVE env had no pandas. create_env never runs again for an env that
        exists, so this step is how it arrives."""
        self.setup()
        os.remove(os.path.join(self.prefix, "conda-meta", "pandas-2.3.2-py311h0000000_0.json"))
        _s, calls = self.setup()
        installs = [c for c in calls if c.startswith("install")]
        self.assertTrue(installs[-1].endswith("-c conda-forge -c bioconda pandas"), installs)
        self.assertNotIn("pandas", [c for c in calls if c.startswith("create")][0])

    def test_check_installs_nothing(self):
        _s, calls = self.setup("--check")
        self.assertEqual(calls, [])
        self.assertEqual(_lines(self.ensure_log), [], "--check ran ensure_dotnet8.sh")

    def test_dotnet8_is_provisioned_and_reported(self):
        s, _calls = self.setup()
        self.assertEqual(_lines(self.ensure_log), ["call"])
        self.assertEqual(s["dotnet8"]["root"], self.dotnet_root)
        self.assertIn(self.dotnet_root, s["dotnet8"]["note"])
        self.assertFalse(any("ensure_dotnet8.sh could not" in n for n in s["notes"]), s["notes"])

    def test_a_failed_dotnet8_is_a_note_not_a_failure(self):
        s, _calls = self.setup(FAKE_ENSURE_FAIL="1")      # setup() asserts exit 0
        self.assertEqual(_lines(self.ensure_log), ["call"])
        self.assertEqual(s["dotnet8"]["root"], "")
        self.assertIn("could not provide .NET 8", s["dotnet8"]["note"])
        self.assertTrue(any("ensure_dotnet8.sh could not provide .NET 8" in n
                            and "DIA-NN cannot read .raw" in n for n in s["notes"]), s["notes"])
        self.assertIs(s["ready_for"]["de"], True, "a .NET failure took something else down")

    def test_a_second_run_reuses_it(self):
        self.setup()
        s, _calls = self.setup()
        self.assertEqual(_lines(self.ensure_log), ["call", "call"])
        self.assertEqual(s["dotnet8"]["root"], self.dotnet_root)

    def test_a_failed_parser_install_costs_nothing_else(self):
        s, _calls = self.setup(FAKE_MM_INSTALL_FAIL="1")
        self.assertIs(s["ready_for"]["de"], True, "R went down with the parser")
        self.assertIs(s["ready_for"]["thermo_raw"], False)
        self.assertTrue(any("Could not install thermorawfileparser pythonnet pandas" in n
                            for n in s["notes"]), s["notes"])


class ActivateExportsTheDotnetRoot(setup_t.SetupCheckHarness):
    def _source(self, **env):
        self.check()
        full = {"PATH": self.sys, "HOME": self.d}
        full.update(env)
        p = subprocess.run([self.bash, "-c", f'. "{self.home}/activate.sh"; '
                            'printf "%s\\n%s\\n" "${DOTNET_ROOT:-}" "$PATH"'],
                           capture_output=True, text=True, env=full)
        self.assertEqual(p.returncode, 0, p.stderr)
        return p.stdout.splitlines()

    def test_ensure_dotnet8s_install_is_exported(self):
        root = fake_dotnet_root(os.path.join(self.d, ".proteomics-pipeline", "dotnet8"))
        dotnet_root, path = self._source()
        self.assertEqual(dotnet_root, root)
        self.assertTrue(path.startswith(root + os.pathsep), path)

    def test_a_dotnet_root_already_set_is_kept(self):
        fake_dotnet_root(os.path.join(self.d, ".proteomics-pipeline", "dotnet8"))
        self.assertEqual(self._source(DOTNET_ROOT="/opt/mine")[0], "/opt/mine")

    def test_nothing_is_exported_before_it_exists(self):
        self.assertEqual(self._source()[0], "")


if __name__ == "__main__":
    unittest.main(verbosity=2)
