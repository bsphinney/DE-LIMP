#!/usr/bin/env python3
"""
setup.json's `diann.ready` is what the orchestrator gates a DIA run on, and `diann.note` is
the instruction a user follows when it is false. Both are decided by one chain in setup.sh.

The chain used to ask for /quobyte/proteomics-grp AND apptainer. The apptainer half was
wrong and was removed: HIVE keeps DIA-NN as NATIVE builds (build_<nnn>/diann-<version>/
diann-linux) and only the old 2.3.0 is a .sif, so reusing them needs no container runtime.

But removing it left `$QUOBYTE` alone, and /quobyte is a network filesystem that can be
mounted anywhere -- including a Mac, where DIA-NN has no native build at all. Such a host
reported `diann.ready: true` with a note telling it to reuse binaries it cannot execute,
instead of the Docker Desktop instructions the macOS branch gives. The test is the OS.

`bash setup.sh --check` installs nothing; PP_HOME and QUOBYTE_DIR are pointed at temp
directories and PATH holds only a handful of linked coreutils, so no conda, Rscript, sage,
docker or apptainer on the machine running this can change the answer.
"""
import json
import os
import platform
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
SETUP = os.path.join(SCRIPTS, "setup.sh")

COREUTILS = ("cat", "mkdir", "sed", "tee", "tr", "uname")
IS_MAC = platform.system() == "Darwin"
IS_LINUX = platform.system() == "Linux"


class SetupCheckHarness(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.sys = os.path.join(self.d, "sysbin")
        self.home = os.path.join(self.d, "pp")
        self.quobyte = os.path.join(self.d, "quobyte-proteomics-grp")
        for p in (self.sys, self.quobyte):
            os.makedirs(p)
        for tool in COREUTILS:
            real = shutil.which(tool)
            if not real:
                self.skipTest(f"{tool} is not installed; setup.sh needs it")
            os.symlink(real, os.path.join(self.sys, tool))
        self.bash = shutil.which("bash")

    def tearDown(self):
        self._tmp.cleanup()

    def check(self, quobyte_dir=None, extra_env=None):
        env = {"PATH": self.sys, "HOME": self.d, "PP_HOME": self.home,
               "QUOBYTE_DIR": self.quobyte if quobyte_dir is None else quobyte_dir}
        env.update(extra_env or {})
        r = subprocess.run([self.bash, SETUP, "--check"], capture_output=True, text=True,
                           env=env, timeout=300)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(self.home, "setup.json")) as fh:
            return json.load(fh)


class QuobyteDirOverrideTests(SetupCheckHarness):
    def test_the_shared_folder_is_where_QUOBYTE_DIR_says(self):
        """Hard-coding the path meant this chain could only ever be tested on HIVE itself."""
        self.assertIs(self.check()["uc_davis_hive"], True)

    def test_a_missing_shared_folder_is_not_hive(self):
        s = self.check(quobyte_dir=os.path.join(self.d, "nope"))
        self.assertIs(s["uc_davis_hive"], False)


class DiannReadinessTests(SetupCheckHarness):
    @unittest.skipUnless(IS_MAC, "the macOS branch")
    def test_a_mac_with_a_quobyte_mount_is_not_diann_ready(self):
        """No docker either: DIA-NN cannot run on this host by any route, and the note has
        to be the one that says how to fix that."""
        s = self.check()
        self.assertIs(s["uc_davis_hive"], True)
        self.assertIs(s["diann"]["ready"], False)
        self.assertIn("macOS", s["diann"]["note"])
        self.assertNotIn("HIVE", s["diann"]["note"])
        self.assertIs(s["ready_for"]["dia"], False)

    @unittest.skipUnless(IS_MAC, "the macOS branch")
    def test_a_mac_without_a_quobyte_mount_says_the_same_thing(self):
        s = self.check(quobyte_dir=os.path.join(self.d, "nope"))
        self.assertIs(s["diann"]["ready"], False)
        self.assertIn("macOS", s["diann"]["note"])

    @unittest.skipUnless(IS_LINUX, "the HIVE branch")
    def test_hive_is_diann_ready_without_apptainer(self):
        """The half of the old condition that WAS wrong: the Core's builds are native
        binaries, so a HIVE login shell with no apptainer on PATH is still ready."""
        s = self.check()
        self.assertIs(s["has_apptainer"], False)
        self.assertIs(s["diann"]["ready"], True)
        self.assertIn("HIVE", s["diann"]["note"])
        self.assertIn(self.quobyte, s["diann"]["note"])

    @unittest.skipUnless(IS_LINUX, "the plain-Linux branch")
    def test_linux_without_the_shared_folder_downloads_instead(self):
        s = self.check(quobyte_dir=os.path.join(self.d, "nope"))
        self.assertIs(s["diann"]["ready"], True)
        self.assertIn("acquire_tools.sh downloads", s["diann"]["note"])


if __name__ == "__main__":
    unittest.main()
