#!/usr/bin/env python3
"""
tools.json `versions` must name the build that runs, never the request.

references/reproducibility.md promises that the recorded version "is the one actually
installed -- never the literal string `latest`", and run_search.py copies it verbatim into
search_provenance.json, which fran_deposit.py hands to FRAN as the engine version. Only the
DIA-NN *download* branch kept that promise. Every other branch that resolves `latest` wrote
the request back out:

  * HIVE native builds: `find ... | sort -V | tail -n1` picked a real binary, then
    recorded "latest".
  * HIVE .sif: same.
  * Sage: the cached tarball under sage/latest/ recorded "latest". Measured on HIVE
    2026-09-16: BOTH the shared ~/.proteomics-pipeline/tools/tools.json and the FRAN pilot's
    tools.json say `"sage": "latest"`; the tarball is sage-v0.14.7, and the binary's own
    `--version` says "sage 0.14.6" -- so nothing in the record could tell anyone which it was.

Everything here runs against a MOCKED layout: a fake HIVE build directory, a fake tools root,
and a stub `curl` that serves canned GitHub listings. PATH holds ONLY that stub directory and
a directory of links to the few coreutils the script calls -- never /usr/bin itself -- so a
machine's own `sage` (SageMath installs /usr/bin/sage on Linux), `docker`, `apptainer` or
`pip` cannot change the answer, and acquire_alphadia's `pip install` cannot fire on a host that
has HIVE's alphadia.sif but no apptainer. No network, no real engine.
"""
import io
import json
import os
import shutil
import stat
import subprocess
import tarfile
import tempfile
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
ACQUIRE = os.path.join(SCRIPTS, "acquire_tools.sh")

DIANN_RELEASES = """[
  {"assets": [
    {"browser_download_url": "https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.6.1-Academia-Linux.zip"},
    {"browser_download_url": "https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.6.1-Academia.msi"},
    {"browser_download_url": "https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.7.0-Academia-Linux.zip"},
    {"browser_download_url": "https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.7.0-Academia.msi"}
  ]}
]
"""

# A curl that knows three URLs and fails (exit 22, like `curl -f`) on everything else, so
# FragPipe / Sage release lookups simply resolve nothing instead of touching the network.
STUB_CURL = r"""#!/bin/sh
url=""; out=""
while [ $# -gt 0 ]; do
  case "$1" in
    -o) out="$2"; shift 2 ;;
    http*) url="$1"; shift ;;
    *) shift ;;
  esac
done
case "$url" in
  *api.github.com/repos/vdemichev/DiaNN/releases*) cat "$FIXTURES/diann_releases.json" ;;
  *DIA-NN-*-Academia-Linux.zip) cp "$FIXTURES/diann.zip" "$out" ;;
  *) exit 22 ;;
esac
"""


# What acquire_tools.sh and diann_release.sh call, plus what STUB_CURL itself needs. Linked
# one by one into the harness PATH; anything else (sage, docker, apptainer, pip, java,
# alphadia) is absent unless a test stubs it.
COREUTILS = ("basename", "cat", "chmod", "cp", "cut", "dirname", "find", "grep", "gzip", "head",
             "ls", "mkdir", "rm", "sed", "sort", "tail", "tar", "tee", "tr", "uname", "unzip",
             "xargs")


def _exe(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


class AcquireHarness(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.bin = os.path.join(self.d, "bin")
        self.fix = os.path.join(self.d, "fixtures")
        self.root = os.path.join(self.d, "tools")
        self.hive = os.path.join(self.d, "dia-nn")
        self.sys = os.path.join(self.d, "sysbin")
        for p in (self.bin, self.fix, self.root, self.hive, self.sys):
            os.makedirs(p)
        for tool in COREUTILS:
            real = shutil.which(tool)
            if not real:
                self.skipTest(f"{tool} is not installed; acquire_tools.sh needs it")
            os.symlink(real, os.path.join(self.sys, tool))
        self.bash = shutil.which("bash")
        _exe(os.path.join(self.bin, "curl"), STUB_CURL)
        with open(os.path.join(self.fix, "diann_releases.json"), "w") as fh:
            fh.write(DIANN_RELEASES)
        with zipfile.ZipFile(os.path.join(self.fix, "diann.zip"), "w") as z:
            info = zipfile.ZipInfo("diann-2.7.0/diann-linux")
            info.external_attr = 0o755 << 16
            z.writestr(info, "#!/bin/sh\necho 'DIA-NN 2.7.0 Academia'\n")

    def tearDown(self):
        self._tmp.cleanup()

    def hive_build(self, rel):
        """A native build at <hive>/<rel>/diann-linux, e.g. build_261/diann-2.6.1."""
        p = os.path.join(self.hive, rel)
        os.makedirs(p, exist_ok=True)
        _exe(os.path.join(p, "diann-linux"), "#!/bin/sh\nexit 0\n")
        return os.path.join(p, "diann-linux")

    def acquire(self, platform_class, pin_engine="", pin_version="", extra_env=None):
        # PATH is the stub dir plus the linked coreutils only (see COREUTILS): no real curl,
        # sage, docker, apptainer or pip on the machine running the tests can leak in.
        env = {"PATH": f"{self.bin}:{self.sys}", "HOME": self.d, "FIXTURES": self.fix,
               "DIANN_HIVE_DIR": self.hive, "PIN_ENGINE": pin_engine,
               "PIN_VERSION": pin_version}
        env.update(extra_env or {})
        r = subprocess.run([self.bash, ACQUIRE, platform_class, self.root],
                           capture_output=True, text=True, env=env, timeout=120)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(self.root, "tools.json")) as fh:
            return json.load(fh)


class HiveNativeBuildTests(AcquireHarness):
    def test_latest_records_the_resolved_native_build_not_the_word_latest(self):
        """The HIVE layout as listed on 2026-09-16: 2.5.1, 2.6.0 and 2.6.1 native builds."""
        for rel in ("build_251/diann-2.5.1", "build_260/diann-2.6.0", "build_261/diann-2.6.1"):
            self.hive_build(rel)
        t = self.acquire("hpc")
        self.assertTrue(t["diann"].endswith("build_261/diann-2.6.1/diann-linux"), t["diann"])
        self.assertEqual(t["versions"]["diann"], "2.6.1")

    def test_latest_is_the_highest_VERSION_not_the_last_path(self):
        """`sort -V` over whole paths orders by the directory names in front of the version.
        A build kept outside the build_<nnn> convention then wins or loses on its folder name:
        here the old 2.5.1 in `z_archive/` sorts after `build_261/` and would be chosen."""
        self.hive_build("build_261/diann-2.6.1")
        self.hive_build("z_archive/diann-2.5.1")
        t = self.acquire("hpc")
        self.assertEqual(t["versions"]["diann"], "2.6.1")
        self.assertIn("diann-2.6.1", t["diann"])

    def test_a_pinned_native_build_still_records_its_pin(self):
        self.hive_build("build_260/diann-2.6.0")
        self.hive_build("build_261/diann-2.6.1")
        t = self.acquire("hpc", "diann", "2.6.0")
        self.assertEqual(t["versions"]["diann"], "2.6.0")
        self.assertIn("diann-2.6.0", t["diann"])

    def test_latest_with_only_a_sif_records_the_sif_version(self):
        with open(os.path.join(self.hive, "diann_2.3.0.sif"), "w") as fh:
            fh.write("not really an image")
        _exe(os.path.join(self.bin, "apptainer"), "#!/bin/sh\nexit 0\n")
        t = self.acquire("hpc")
        self.assertIn("diann_2.3.0.sif", t["diann"])
        self.assertEqual(t["versions"]["diann"], "2.3.0")


class DownloadBranchTests(AcquireHarness):
    def test_latest_download_records_the_downloaded_build(self):
        """Already correct before this change (the download branch resolves by asset
        filename); kept as a guard so the refactor cannot regress it."""
        t = self.acquire("linux")
        self.assertEqual(t["versions"]["diann"], "2.7.0")
        self.assertTrue(t["diann"].endswith("diann/2.7.0/diann-2.7.0/diann-linux"), t["diann"])

    def test_a_hive_pin_that_is_absent_downloads_and_records_that_pin(self):
        """The FRAN pilot's path: pin 2.7.0 on HIVE, where only 2.5.1-2.6.1 exist."""
        self.hive_build("build_261/diann-2.6.1")
        t = self.acquire("hpc", "diann", "2.7.0")
        self.assertEqual(t["versions"]["diann"], "2.7.0")
        self.assertNotIn("2.6.1", t["diann"])

    def test_nothing_resolvable_records_no_version_rather_than_latest(self):
        """No HIVE build and no reachable release: there is no command, so there is no
        version either -- an empty string, not the request echoed back."""
        os.remove(os.path.join(self.fix, "diann_releases.json"))
        t = self.acquire("hpc")
        self.assertIsNone(t["diann"])
        self.assertNotEqual(t["versions"]["diann"], "latest")
        self.assertEqual(t["versions"]["diann"], "")


class DockerImageTests(AcquireHarness):
    def test_mac_docker_image_records_the_version_in_its_tag(self):
        """build_diann_docker.sh tags the image with the RESOLVED version, so the tag is the
        record; `latest` must not be written over it."""
        _exe(os.path.join(self.bin, "docker"), "#!/bin/sh\nexit 0\n")
        t = self.acquire("mac", extra_env={"DIANN_DOCKER_IMAGE": "proteomics-pipeline/diann:2.6.1"})
        self.assertIn("proteomics-pipeline/diann:2.6.1", t["diann"])
        self.assertEqual(t["versions"]["diann"], "2.6.1")

    def test_an_untagged_image_records_no_version_and_says_why(self):
        _exe(os.path.join(self.bin, "docker"), "#!/bin/sh\nexit 0\n")
        t = self.acquire("mac", extra_env={"DIANN_DOCKER_IMAGE": "mylab/diann"})
        self.assertEqual(t["versions"]["diann"], "")
        self.assertTrue(any("mylab/diann" in n and "version" in n for n in t["notes"]), t["notes"])


class SageTests(AcquireHarness):
    def _cached_sage(self, cache_name, tag="v0.14.7"):
        """What acquire_sage leaves behind: the binary, stripped out of its top-level folder,
        next to the tarball it came from. The tarball keeps the release folder name."""
        d = os.path.join(self.root, "sage", cache_name)
        os.makedirs(d)
        # The real v0.14.7 binary reports 0.14.6 (measured on HIVE). The stub does the same,
        # so a version read from `--version` would visibly be the wrong answer here.
        _exe(os.path.join(d, "sage"), "#!/bin/sh\necho 'sage 0.14.6'\n")
        top = f"sage-{tag}-x86_64-unknown-linux-gnu"
        buf = io.BytesIO(b"#!/bin/sh\n")
        with tarfile.open(os.path.join(d, "sage.tar.gz"), "w:gz") as tf:
            ti = tarfile.TarInfo(top)
            ti.type = tarfile.DIRTYPE
            tf.addfile(ti)
            ti = tarfile.TarInfo(f"{top}/sage")
            ti.size = len(buf.getvalue())
            tf.addfile(ti, buf)
        return d

    def test_a_cached_latest_sage_records_its_release_not_latest(self):
        self._cached_sage("latest")
        t = self.acquire("hpc")
        self.assertTrue(t["sage"].endswith("sage/latest/sage"), t["sage"])
        self.assertEqual(t["versions"]["sage"], "0.14.7")

    def test_sage_version_comes_from_the_release_not_the_binarys_self_report(self):
        """`sage --version` of the v0.14.7 release prints 0.14.6. The release tag is what
        resolve_defaults.py pins (0.14.7) and what anyone can download again."""
        self._cached_sage("latest")
        t = self.acquire("linux")
        self.assertEqual(t["versions"]["sage"], "0.14.7")     # not 0.14.6, and not "" either

    def test_a_sage_on_PATH_is_not_a_release(self):
        """A `sage` on PATH records "env": a source, not a build. run_search.py must then
        record no version -- see test_engine_version_provenance.py."""
        _exe(os.path.join(self.bin, "sage"), "#!/bin/sh\necho 'sage 0.14.6'\n")
        t = self.acquire("linux")
        self.assertEqual(t["sage"], os.path.join(self.bin, "sage"))
        self.assertEqual(t["versions"]["sage"], "env")


class RadiantTests(AcquireHarness):
    """Radiant is only acquired on request (ACQUIRE_RADIANT=1 or PIN_ENGINE=radiant)."""

    def test_an_unpinned_radiant_image_records_no_version_not_latest(self):
        """`seerbio/radiant-fulcrum:latest` names no release, and resolving what :latest
        points at would mean pulling ~3 GB. So nothing is recorded, and a note says how to
        get a recorded release (pin one)."""
        _exe(os.path.join(self.bin, "docker"), "#!/bin/sh\nexit 0\n")
        t = self.acquire("linux", extra_env={"ACQUIRE_RADIANT": "1"})
        self.assertEqual(t["radiant_image"], "seerbio/radiant-fulcrum:latest")
        self.assertEqual(t["versions"]["radiant"], "")
        self.assertTrue(any("PIN_ENGINE=radiant" in n for n in t["notes"]), t["notes"])

    def test_a_pinned_radiant_image_records_its_tag(self):
        _exe(os.path.join(self.bin, "docker"), "#!/bin/sh\nexit 0\n")
        t = self.acquire("linux", "radiant", "2.3.3")
        self.assertEqual(t["radiant_image"], "seerbio/radiant-fulcrum:2.3.3")
        self.assertEqual(t["versions"]["radiant"], "2.3.3")

    def test_no_container_runtime_records_no_version(self):
        t = self.acquire("linux", "radiant", "2.3.3")
        self.assertIsNone(t["radiant"])
        self.assertEqual(t["versions"]["radiant"], "")


if __name__ == "__main__":
    unittest.main()
