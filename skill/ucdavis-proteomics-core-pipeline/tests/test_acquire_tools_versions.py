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
  # FragPipe only answers for the tests that lay the fixtures down; everything else still
  # resolves nothing, exactly as before.
  *api.github.com/repos/Nesvilab/FragPipe/releases*)
     [ -f "$FIXTURES/fragpipe_releases.json" ] || exit 22
     cat "$FIXTURES/fragpipe_releases.json" ;;
  *FragPipe-*.zip)
     [ -f "$FIXTURES/fragpipe.zip" ] || exit 22
     cp "$FIXTURES/fragpipe.zip" "$out" ;;
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
        # stand-ins for /quobyte/proteomics-grp/apptainers and /quobyte/proteomics-grp/radiant
        self.radiant_dirs = (os.path.join(self.d, "apptainers"), os.path.join(self.d, "radiant"))
        for p in (self.bin, self.fix, self.root, self.hive, self.sys, *self.radiant_dirs):
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
        # RADIANT_HIVE_DIRS points the HIVE Radiant .sif search at empty mocked folders, so
        # the real /quobyte/proteomics-grp/apptainers on HIVE cannot answer for a test.
        env = {"PATH": f"{self.bin}:{self.sys}", "HOME": self.d, "FIXTURES": self.fix,
               "DIANN_HIVE_DIR": self.hive, "PIN_ENGINE": pin_engine,
               "PIN_VERSION": pin_version,
               "RADIANT_HIVE_DIRS": f"{self.radiant_dirs[0]}:{self.radiant_dirs[1]}"}
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


class SageReleaseFallbackTests(AcquireHarness):
    """The release recorded for Sage is the one the cached TARBALL names, or nothing.

    acquire_sage caches under sage/<REQUEST>/ and the binary cannot be asked (the v0.14.7
    release binary prints 0.14.6). When the tarball is gone -- copied installs, a cache
    someone cleaned -- there is nothing left that names the release. The fallback used to
    record the REQUEST instead, i.e. PIN_VERSION, i.e. the workflow manifest's own pin: the
    build record then said what the manifest asked for, and run_search.py compared the pin
    with itself (the cache path sage/<pin>/sage names no version it will read) and wrote
    `mismatch: false` -- a claim that the two were checked and agree."""

    _cached_sage = SageTests._cached_sage          # the same layout, tarball and all

    def _binary_only(self, cache_name):
        d = os.path.join(self.root, "sage", cache_name)
        os.makedirs(d)
        _exe(os.path.join(d, "sage"), "#!/bin/sh\necho 'sage 0.14.6'\n")
        return d

    def test_a_pinned_sage_with_no_tarball_records_no_version_not_the_pin(self):
        d = self._binary_only("0.14.7")
        t = self.acquire("linux", "sage", "0.14.7")
        self.assertTrue(t["sage"].endswith("sage/0.14.7/sage"), t["sage"])
        self.assertEqual(t["versions"]["sage"], "")
        self.assertTrue(any("cannot tell which release" in n and d in n for n in t["notes"]),
                        t["notes"])

    def test_an_unpinned_sage_with_no_tarball_records_no_version_either(self):
        self._binary_only("latest")
        t = self.acquire("linux")
        self.assertEqual(t["versions"]["sage"], "")

    def test_a_pinned_sage_still_records_the_release_the_tarball_names(self):
        """The evidence, when it is there, is still read -- and when it disagrees with the
        pin, the CACHE wins and the disagreement is a note."""
        self._cached_sage("0.14.7", tag="v0.14.6")
        t = self.acquire("linux", "sage", "0.14.7")
        self.assertEqual(t["versions"]["sage"], "0.14.6")
        self.assertTrue(any("pinned '0.14.7'" in n and "0.14.6" in n for n in t["notes"]),
                        t["notes"])


class FragPipeVersionTests(AcquireHarness):
    """tools.json must carry a `fragpipe` version. Without the key, run_search.py has nothing
    to read and EVERY FragPipe search records `version: null` for ever -- the engine is
    invisible in the deposited record, not merely unresolved this once."""

    RELEASES = """[
  {"assets": [
    {"browser_download_url": "https://github.com/Nesvilab/FragPipe/releases/download/24.0/FragPipe-24.0.zip"},
    {"browser_download_url": "https://github.com/Nesvilab/FragPipe/releases/download/24.0/FragPipe-jre-24.0.zip"}
  ]}
]
"""

    def fragpipe_release(self, top):
        """The release listing plus a zip that unpacks to <top>/bin/fragpipe."""
        with open(os.path.join(self.fix, "fragpipe_releases.json"), "w") as fh:
            fh.write(self.RELEASES)
        with zipfile.ZipFile(os.path.join(self.fix, "fragpipe.zip"), "w") as z:
            info = zipfile.ZipInfo(f"{top}/bin/fragpipe")
            info.external_attr = 0o755 << 16
            z.writestr(info, "#!/bin/sh\nexit 0\n")

    def test_a_downloaded_fragpipe_records_the_release_it_unpacked(self):
        self.fragpipe_release("fragpipe-24.0")
        t = self.acquire("linux")
        self.assertTrue(t["fragpipe"].endswith("fragpipe-24.0/bin/fragpipe"), t["fragpipe"])
        self.assertEqual(t["versions"]["fragpipe"], "24.0")

    def test_a_folder_that_does_not_name_the_release_falls_back_to_the_asset(self):
        """Some FragPipe zips unpack into a plain `fragpipe/`. The asset that was downloaded
        (FragPipe-24.0.zip) still names the release; the cache folder -- named after the
        REQUEST, `latest` -- never does."""
        self.fragpipe_release("fragpipe")
        t = self.acquire("linux")
        self.assertTrue(t["fragpipe"].endswith("fragpipe/bin/fragpipe"), t["fragpipe"])
        self.assertEqual(t["versions"]["fragpipe"], "24.0")

    def test_the_request_keyed_cache_folder_is_never_read_as_the_release(self):
        """A pinned install caches under fragpipe/24.0/. That folder is named after what was
        ASKED for; if the zip inside it is some other release, the path must not claim 24.0."""
        self.fragpipe_release("fragpipe")
        with zipfile.ZipFile(os.path.join(self.fix, "fragpipe.zip"), "w") as z:
            info = zipfile.ZipInfo("fragpipe/bin/fragpipe")
            info.external_attr = 0o755 << 16
            z.writestr(info, "#!/bin/sh\nexit 0\n")
        os.remove(os.path.join(self.fix, "fragpipe_releases.json"))   # nothing resolves
        os.makedirs(os.path.join(self.root, "fragpipe", "24.0", "bin"))
        fp = os.path.join(self.root, "fragpipe", "24.0", "bin", "fragpipe")
        _exe(fp, "#!/bin/sh\nexit 0\n")
        t = self.acquire("linux", "fragpipe", "24.0")
        self.assertEqual(t["fragpipe"], fp)
        self.assertEqual(t["versions"]["fragpipe"], "")
        self.assertTrue(any("cannot tell which release" in n for n in t["notes"]), t["notes"])

    def test_no_fragpipe_at_all_records_an_empty_version_not_a_missing_key(self):
        t = self.acquire("linux")
        self.assertIsNone(t["fragpipe"])
        self.assertEqual(t["versions"]["fragpipe"], "")


class AlphaDiaVersionTests(AcquireHarness):
    """AlphaDIA is a pip package, so pip's own metadata is the record. Without a `versions`
    key it was the same permanent `version: null` as FragPipe."""

    def stub_alphadia(self, pip_version="1.10.0"):
        _exe(os.path.join(self.bin, "alphadia"), "#!/bin/sh\nexit 0\n")
        _exe(os.path.join(self.bin, "pip"),
             "#!/bin/sh\n"
             "[ \"$1\" = show ] && [ \"$2\" = alphadia ] || exit 1\n"
             f"echo 'Name: alphadia'\necho 'Version: {pip_version}'\n")

    def test_an_installed_alphadia_records_what_pip_reports(self):
        self.stub_alphadia()
        t = self.acquire("linux")
        self.assertEqual(t["alphadia"], os.path.join(self.bin, "alphadia"))
        self.assertEqual(t["versions"]["alphadia"], "1.10.0")

    def test_a_prerelease_is_recorded_verbatim_not_rounded_to_a_release(self):
        """run_search.py refuses it as `version` (it is not shaped like a release), which is
        the point: 1.10.0rc1 is not 1.10.0, and neither is a lie about the other."""
        self.stub_alphadia("1.10.0rc1")
        t = self.acquire("linux")
        self.assertEqual(t["versions"]["alphadia"], "1.10.0rc1")

    def test_an_alphadia_pip_cannot_account_for_records_no_version_and_says_so(self):
        _exe(os.path.join(self.bin, "alphadia"), "#!/bin/sh\nexit 0\n")
        _exe(os.path.join(self.bin, "pip"), "#!/bin/sh\nexit 1\n")
        t = self.acquire("linux")
        self.assertEqual(t["versions"]["alphadia"], "")
        self.assertTrue(any("pip reports no version" in n for n in t["notes"]), t["notes"])

    def test_no_alphadia_records_an_empty_version_not_a_missing_key(self):
        t = self.acquire("linux")
        self.assertIsNone(t["alphadia"])
        self.assertEqual(t["versions"]["alphadia"], "")


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

    def test_a_pinned_image_says_the_tag_is_all_the_evidence_there_is(self):
        """The version here is read off the image REFERENCE that was recorded and will be
        run, not off $ver -- but nothing has pulled that image, so the tag is the only thing
        that says what is inside it. Two sources that are one source: say so in the note."""
        _exe(os.path.join(self.bin, "docker"), "#!/bin/sh\nexit 0\n")
        t = self.acquire("linux", "radiant", "2.3.3")
        self.assertTrue(any("only evidence of what is inside" in n for n in t["notes"]),
                        t["notes"])

    def test_an_apptainer_only_host_records_the_tag_the_same_way(self):
        _exe(os.path.join(self.bin, "apptainer"), "#!/bin/sh\nexit 0\n")
        t = self.acquire("linux", "radiant", "2.3.3")
        self.assertEqual(t["radiant_image"], "docker://seerbio/radiant-fulcrum:2.3.3")
        self.assertEqual(t["versions"]["radiant"], "2.3.3")
        t = self.acquire("linux", extra_env={"ACQUIRE_RADIANT": "1"})
        self.assertEqual(t["versions"]["radiant"], "")

    def test_no_container_runtime_records_no_version(self):
        t = self.acquire("linux", "radiant", "2.3.3")
        self.assertIsNone(t["radiant"])
        self.assertEqual(t["versions"]["radiant"], "")


class HiveRadiantSifTests(AcquireHarness):
    """On HIVE Radiant is a prebuilt .sif, found in /quobyte/proteomics-grp/apptainers or
    /quobyte/proteomics-grp/radiant (listed 2026-09-16: apptainers/radiant-fulcrum-2.3.3.sif).
    The release recorded is the one the .sif's NAME carries; RADIANT_HIVE_DIRS mocks the two
    folders here, as DIANN_HIVE_DIR does for DIA-NN."""

    def setUp(self):
        super().setUp()
        _exe(os.path.join(self.bin, "apptainer"), "#!/bin/sh\nexit 0\n")

    def sif(self, name, where=0):
        p = os.path.join(self.radiant_dirs[where], name)
        with open(p, "w") as fh:
            fh.write("not really an image")
        return p

    def test_an_unpinned_hive_sif_records_the_release_its_name_carries(self):
        """HIVE's real layout. origin/main recorded "latest" here (srun 23512217)."""
        p = self.sif("radiant-fulcrum-2.3.3.sif")
        t = self.acquire("hpc", extra_env={"ACQUIRE_RADIANT": "1"})
        self.assertEqual(t["radiant_image"], p)
        self.assertEqual(t["versions"]["radiant"], "2.3.3")

    def test_unpinned_takes_the_highest_release_by_version(self):
        self.sif("radiant-fulcrum-2.3.2.sif")
        p = self.sif("radiant-fulcrum-2.3.10.sif")
        t = self.acquire("hpc", extra_env={"ACQUIRE_RADIANT": "1"})
        self.assertEqual(t["radiant_image"], p)
        self.assertEqual(t["versions"]["radiant"], "2.3.10")

    def test_a_pinned_hive_sif_records_that_release(self):
        self.sif("radiant-fulcrum-2.3.2.sif")
        p = self.sif("radiant-fulcrum-2.3.3.sif")
        t = self.acquire("hpc", "radiant", "2.3.3")
        self.assertEqual(t["radiant_image"], p)
        self.assertEqual(t["versions"]["radiant"], "2.3.3")

    def test_the_second_folder_is_searched_too(self):
        p = self.sif("radiant-fulcrum-2.3.3.sif", where=1)
        t = self.acquire("hpc", extra_env={"ACQUIRE_RADIANT": "1"})
        self.assertEqual(t["radiant_image"], p)
        self.assertEqual(t["versions"]["radiant"], "2.3.3")

    def test_a_pinned_sif_whose_name_ends_otherwise_falls_back_to_the_pin_it_matched(self):
        """The name must end `-<ver>.sif` to be read; a pinned search only matches names that
        contain the pin, so the pin is then what the name carries."""
        p = self.sif("radiant-2.3.3-fulcrum.sif")
        t = self.acquire("hpc", "radiant", "2.3.3")
        self.assertEqual(t["radiant_image"], p)
        self.assertEqual(t["versions"]["radiant"], "2.3.3")

    def test_an_unpinned_sif_that_names_no_release_records_none_and_says_so(self):
        p = self.sif("radiant-fulcrum.sif")
        t = self.acquire("hpc", extra_env={"ACQUIRE_RADIANT": "1"})
        self.assertEqual(t["radiant_image"], p)
        self.assertEqual(t["versions"]["radiant"], "")
        self.assertTrue(any("names no release" in n for n in t["notes"]), t["notes"])

    def test_no_hive_sif_records_no_version_and_names_the_folders_searched(self):
        t = self.acquire("hpc", extra_env={"ACQUIRE_RADIANT": "1"})
        self.assertIsNone(t["radiant"])
        self.assertEqual(t["versions"]["radiant"], "")
        self.assertTrue(any(self.radiant_dirs[0] in n for n in t["notes"]), t["notes"])

    def test_an_empty_RADIANT_HIVE_DIRS_searches_nowhere_rather_than_quobyte(self):
        """`${RADIANT_HIVE_DIRS:-<the /quobyte paths>}` made an explicitly EMPTY value mean
        "use HIVE's own folders" -- the one value a test or another site would set to mean
        the opposite. It also has to survive `set -u` on bash 3.2, where expanding an empty
        array aborts the script."""
        t = self.acquire("hpc", "radiant", "2.3.3", extra_env={"RADIANT_HIVE_DIRS": ""})
        self.assertEqual(t["versions"]["radiant"], "")
        self.assertFalse([n for n in t["notes"] if "/quobyte" in n], t["notes"])
        self.assertTrue(any("RADIANT_HIVE_DIRS is empty" in n for n in t["notes"]), t["notes"])

    def test_empty_folders_in_the_list_are_skipped_not_searched(self):
        p = self.sif("radiant-fulcrum-2.3.3.sif", where=1)
        dirs = f"::{self.radiant_dirs[0]}::{self.radiant_dirs[1]}:"
        t = self.acquire("hpc", extra_env={"ACQUIRE_RADIANT": "1", "RADIANT_HIVE_DIRS": dirs})
        self.assertEqual(t["radiant_image"], p)
        self.assertEqual(t["versions"]["radiant"], "2.3.3")


if __name__ == "__main__":
    unittest.main()
