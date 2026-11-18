#!/usr/bin/env python3
"""
tools.json `fragpipe_tools`: whether FragPipe can search here, from what is actually installed
(2.15, lane F2; C0 2026-11-04, Rulings 86/90).

What 2.14.2 got wrong, measured on HIVE by C0:
  * tools.json's `fragpipe` was a home copy whose tools/ has no MSFragger/IonQuant jar; the Core's
    licensed install (/quobyte/proteomics-grp/fragpipe24/fragpipe-24.0) was reached only through
    FRAGPIPE_TOOLS_FOLDER (unset on HIVE) or the shared ~/.config cache.
  * acquire_tools.sh looked for MSFragger*.jar at the top of the tools folder only, case-sensitive;
    the Core's tools_licensed keeps diatracer-2.2.1/diaTracer-2.2.1.jar one level down, capital T,
    and diatracer_parallel.find_jar's `**/diatracer*.jar` glob missed it.
  * Philosopher and the bundled DIA-NN are the RUNNING binary's own (FragPipe ignores the tools
    folder for Philosopher), so they are recorded from that install.
  * FragPipe 24.0 needs Java 11+: Java 1.8 is unavailable, 25 (HIVE's nodes) available.

Fake installs in temporary folders, a fake `java`; acquire_tools.sh in the mocked harness of
test_acquire_tools_versions.py, with python3 added. No network, no real FragPipe.
"""
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest

from job_env import no_hive_marker

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

try:
    import engine_choice as ec  # noqa: E402
except ImportError:      # a fail-first run on the base tree: each test that needs it fails alone
    ec = None
# the environment names engine_choice reads (spelled here so the base tree's tests still run)
CORE_ENV, TOOLS_ENV = "FRAGPIPE_CORE_INSTALL", "FRAGPIPE_TOOLS_FOLDER"
from test_acquire_tools_versions import AcquireHarness, _exe  # noqa: E402

JAVA25 = 'openjdk version "25.0.4.1" 2025-10-21'
JAVA8 = 'openjdk version "1.8.0_392"'
# the licence note acquire_tools.sh has always written (Ruling 56: licence notes unchanged)
LICENCE_NOTE = ("FragPipe {ver}: downloaded, but MSFragger/IonQuant are license-gated and absent. "
                "Run the GUI once to accept the academic license, OR set FRAGPIPE_TOOLS_FOLDER. "
                "Headless search disabled until then.")


def touch(path, text="", mode=0o644):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, mode)
    return path


def fake_install(root, jars=True, nested=False):
    """A FragPipe 24.0 install laid out as the Core's (C0 item 1) under <root>/fragpipe-24.0."""
    inst = os.path.join(root, "fragpipe-24.0")
    touch(os.path.join(inst, "bin", "fragpipe"), "#!/bin/sh\nexit 0\n", 0o755)
    tools = os.path.join(inst, "tools")
    if jars:
        sub = os.path.join(tools, "msfragger-4.4.1") if nested else tools
        touch(os.path.join(sub, "MSFragger-4.4.1.jar"))
        touch(os.path.join(tools, "IonQuant-1.11.20.jar"))
        touch(os.path.join(tools, "diatracer-2.2.1", "diaTracer-2.2.1.jar"))
    touch(os.path.join(tools, "Philosopher", "philosopher-v5.1.3-RC9"), "#!/bin/sh\n", 0o755)
    touch(os.path.join(tools, "diann", "1.8.2_beta_8", "linux", "diann-1.8.1.8"), "#!/bin/sh\n", 0o755)
    touch(os.path.join(tools, "diann", "1.8.2_beta_8", "linux", "libtimsdata.so"))
    return inst


class JarsAndInstalls(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self._env = os.environ.get(CORE_ENV)
        os.environ[CORE_ENV] = os.path.join(self.d, "core", "fragpipe-24.0")

    def tearDown(self):
        if self._env is None:
            os.environ.pop(CORE_ENV, None)
        else:
            os.environ[CORE_ENV] = self._env
        for root, dirs, _ in os.walk(self.d):
            for x in dirs:
                os.chmod(os.path.join(root, x), 0o755)
        self._tmp.cleanup()

    def test_a_versioned_jar_subfolder_is_found_any_case(self):
        """tools_licensed/diatracer-2.2.1/diaTracer-2.2.1.jar (capital T, one level down)."""
        inst = fake_install(os.path.join(self.d, "core"), nested=True)
        tools = os.path.join(inst, "tools")
        self.assertEqual(ec.find_jar(tools, "diatracer"),
                         {"path": os.path.join(tools, "diatracer-2.2.1", "diaTracer-2.2.1.jar"),
                          "version": "2.2.1"})
        self.assertEqual(ec.find_jar(tools, "msfragger")["version"], "4.4.1")

    def test_three_levels_down_is_not_searched_and_the_highest_version_wins(self):
        t = os.path.join(self.d, "t")
        touch(os.path.join(t, "a", "b", "IonQuant-9.9.9.jar"))
        touch(os.path.join(t, "IonQuant-1.9.0.jar"))
        touch(os.path.join(t, "x", "ionquant-1.11.20.jar"))
        self.assertEqual(ec.find_jar(t, "ionquant")["version"], "1.11.20")

    def test_the_core_install_is_licensed_and_names_what_it_runs(self):
        inst = fake_install(os.path.join(self.d, "core"))
        ft = ec.fragpipe_tools(java=ec.java_record(JAVA25.replace('"', "")), pin="24.0")
        self.assertEqual((ft["licensed"], ft["ready"], ft["source"]), (True, True, "core_install"))
        self.assertEqual(ft["binary"], os.path.join(inst, "bin", "fragpipe"))
        self.assertEqual(ft["folder"], os.path.join(inst, "tools"))
        self.assertEqual({k: (v or {}).get("version") for k, v in ft["jars"].items()},
                         {"msfragger": "4.4.1", "ionquant": "1.11.20", "diatracer": "2.2.1"})
        self.assertEqual(ft["philosopher"]["version"], "5.1.3-RC9")
        self.assertEqual(ft["diann"], {"path": os.path.join(inst, "tools", "diann", "1.8.2_beta_8",
                                                            "linux", "diann-1.8.1.8"),
                                       "version": "1.8.2 beta 8"})
        self.assertEqual(ft["java"]["major"], 25)

    @unittest.skipIf(hasattr(os, "geteuid") and os.geteuid() == 0, "root reads a mode-000 folder")
    def test_an_unreadable_install_is_skipped_with_a_reason(self):
        inst = fake_install(os.path.join(self.d, "core"))
        home = fake_install(os.path.join(self.d, "home"), jars=False)
        os.chmod(inst, 0)
        ft = ec.fragpipe_tools(binary=os.path.join(home, "bin", "fragpipe"),
                               java=ec.java_record(JAVA25))
        self.assertIs(ft["licensed"], False)
        self.assertTrue(any("not readable" in s["why"] for s in ft["skipped"]), ft["skipped"])
        self.assertIn("FRAGPIPE_TOOLS_FOLDER", ft["why"])
        # the home copy's own Philosopher is what would run
        self.assertTrue(ft["philosopher"]["path"].startswith(home))

    def test_FRAGPIPE_TOOLS_FOLDER_is_looked_at_first(self):
        mine = os.path.join(self.d, "mine")
        touch(os.path.join(mine, "MSFragger-4.3.jar"))
        touch(os.path.join(mine, "IonQuant-1.10.27.jar"))
        fake_install(os.path.join(self.d, "core"))
        ft = ec.fragpipe_tools(env_folder=mine, java=ec.java_record(JAVA25))
        self.assertEqual((ft["folder"], ft["source"]), (mine, "FRAGPIPE_TOOLS_FOLDER"))

    def test_a_core_install_of_another_release_is_not_substituted_for_the_pin(self):
        fake_install(os.path.join(self.d, "core"))
        b, why = ec.core_fragpipe_binary(pin="23.0")
        self.assertIsNone(b)
        self.assertIn("not substituted", why)

    def test_fragpipe_ready_reads_tools_json_only(self):
        self.assertIsNone(ec.fragpipe_ready({"fragpipe": "/x"})["ready"])
        self.assertIs(ec.fragpipe_ready({"fragpipe_tools": {"ready": False, "why": "w"}})["ready"],
                      False)


class EasyPQPPython(unittest.TestCase):
    """fragpipe_tools.python: the Python FragPipe hands EasyPQP (--config-python) for a spectral-
    library workflow. With a per-job settings folder (F3) FragPipe no longer reads the shared cache's
    bin-python, so tools.json has to name it. Read from files only (pyvenv.cfg, dist-info): nothing
    imported on a login node."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self._old = {k: os.environ.get(k) for k in (CORE_ENV, "FRAGPIPE_PYTHON", "XDG_CONFIG_HOME")}
        os.environ[CORE_ENV] = os.path.join(self.d, "fragpipe24", "fragpipe-24.0")
        os.environ.pop("FRAGPIPE_PYTHON", None)
        os.environ["XDG_CONFIG_HOME"] = os.path.join(self.d, "xdg")

    def tearDown(self):
        for k, v in self._old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        for root, dirs, _ in os.walk(self.d):
            for x in dirs:
                os.chmod(os.path.join(root, x), 0o755)
        self._tmp.cleanup()

    def venv(self, root, target=sys.executable):
        """A venv as the Core's easypqp_venv2 (2026-11-04): bin/python a link to a conda python."""
        touch(os.path.join(root, "pyvenv.cfg"), "home = /x/bin\nversion = 3.10.17\n")
        os.makedirs(os.path.join(root, "lib", "python3.10", "site-packages",
                                 "easypqp-0.1.59.dist-info"))
        os.makedirs(os.path.join(root, "bin"))
        os.symlink(target, os.path.join(root, "bin", "python"))
        return os.path.join(root, "bin", "python")

    def test_the_core_venv_is_found_with_its_versions(self):
        py = self.venv(os.path.join(self.d, "fragpipe24", "easypqp_venv2"))
        rec = ec.easypqp_python("24.0")
        self.assertEqual({k: rec[k] for k in ("path", "version", "easypqp", "source", "runnable", "why")},
                         {"path": py, "version": "3.10.17", "easypqp": "0.1.59",
                          "source": "core_install", "runnable": True, "why": None})

    @unittest.skipIf(hasattr(os, "geteuid") and os.geteuid() == 0, "root reads a mode-000 folder")
    def test_a_python_linked_into_a_home_this_account_cannot_open_is_not_runnable(self):
        home = os.path.join(self.d, "home_of_someone")
        conda_py = touch(os.path.join(home, "envs", "x", "bin", "python"), "#!/bin/sh\n", 0o755)
        self.venv(os.path.join(self.d, "fragpipe24", "easypqp_venv2"), target=conda_py)
        os.chmod(home, 0)
        rec = ec.easypqp_python("24.0")
        self.assertIs(rec["runnable"], False)
        self.assertIn("cannot run", rec["why"])

    def test_FRAGPIPE_PYTHON_then_the_settings_cache_then_the_core(self):
        self.venv(os.path.join(self.d, "fragpipe24", "easypqp_venv2"))
        cache_py = self.venv(os.path.join(self.d, "mine_venv"))
        touch(os.path.join(self.d, "xdg", "FragPipe", "fragpipe", "24.0", "fragpipe-ui.cache"),
              f"fragpipe-config.bin-python={cache_py}\n")
        self.assertEqual(ec.easypqp_python("24.0")["source"], "fragpipe settings cache")
        env_py = self.venv(os.path.join(self.d, "env_venv"))
        os.environ["FRAGPIPE_PYTHON"] = env_py
        self.assertEqual(ec.easypqp_python("24.0")["path"], env_py)

    def test_fragpipe_tools_carries_it_and_none_when_there_is_none(self):
        fake_install(os.path.join(self.d, "fragpipe24"))
        ft = ec.fragpipe_tools(java=ec.java_record(JAVA25), pin="24.0")
        self.assertIn("python", ft)
        self.assertIsNone(ft["python"])
        self.venv(os.path.join(self.d, "fragpipe24", "easypqp_venv2"))
        self.assertEqual(ec.fragpipe_tools(java=ec.java_record(JAVA25), pin="24.0")["python"]["easypqp"],
                         "0.1.59")


class Java(unittest.TestCase):
    def test_java_1_8_is_unavailable_and_25_available(self):
        for line, major, ok in ((JAVA8, 8, False), (JAVA8.replace('"', ""), 8, False),
                                (JAVA25, 25, True), (JAVA25.replace('"', ""), 25, True),
                                ('openjdk version "11.0.2" 2019-01-15', 11, True),
                                ('java version "17-ea"', 17, True)):
            j = ec.java_record(line)
            self.assertEqual((j["major"], j["ok"]), (major, ok), line)
        self.assertIs(ec.java_record(None)["ok"], False)
        self.assertIsNone(ec.java_record("Unable to locate a Java Runtime.")["ok"])

    def test_java_status_runs_the_java_on_path(self):
        with tempfile.TemporaryDirectory() as d:
            java = _write_java(d, JAVA8)
            self.assertEqual(ec.java_status(java)["major"], 8)

    def test_detect_env_marks_fragpipe_unavailable_under_java_8_and_available_under_25(self):
        for line, ok in ((JAVA8, False), (JAVA25, True)):
            with self.subTest(line), tempfile.TemporaryDirectory() as d:
                b = os.path.join(d, "bin")
                os.makedirs(b)
                _write_java(b, line)
                _exe(os.path.join(b, "uname"),
                     '#!/bin/sh\ncase "$1" in -s) echo Linux;; -m) echo x86_64;; esac\n')
                os.symlink(sys.executable, os.path.join(b, "python3"))
                for t in ("head", "tr", "sed", "dirname", "cat"):
                    os.symlink(shutil.which(t), os.path.join(b, t))
                r = subprocess.run([shutil.which("bash"), os.path.join(SCRIPTS, "detect_env.sh")],
                                   capture_output=True, text=True, env={"PATH": b, **no_hive_marker(d)}, timeout=60)
                self.assertEqual(r.returncode, 0, r.stderr)
                fp = json.loads(r.stdout)["engines"]["fragpipe"]
                self.assertIs(fp["available"], ok)
                if not ok:
                    self.assertIn("Java 8 found", fp["note"])

    def test_detect_env_on_a_mac_stays_unavailable(self):
        with tempfile.TemporaryDirectory() as d:
            b = os.path.join(d, "bin")
            os.makedirs(b)
            _write_java(b, JAVA25)
            _exe(os.path.join(b, "uname"),
                 '#!/bin/sh\ncase "$1" in -s) echo Darwin;; -m) echo arm64;; esac\n')
            for t in ("head", "tr", "sed", "dirname", "cat"):
                os.symlink(shutil.which(t), os.path.join(b, t))
            r = subprocess.run([shutil.which("bash"), os.path.join(SCRIPTS, "detect_env.sh")],
                               capture_output=True, text=True, env={"PATH": b, **no_hive_marker(d)}, timeout=60)
            fp = json.loads(r.stdout)["engines"]["fragpipe"]
            self.assertIs(fp["available"], False)
            self.assertIn("NO macOS build", fp["note"])


def _write_java(d, line):
    p = os.path.join(d, "java")
    _exe(p, f"#!/bin/sh\necho '{line}' >&2\n")
    return p


class AcquireFragPipe(AcquireHarness):
    """acquire_tools.sh with python3 and a Java on PATH, and a mocked Core install."""

    def setUp(self):
        super().setUp()
        self.core = os.path.join(self.d, "core")
        os.symlink(sys.executable, os.path.join(self.bin, "python3"))
        _write_java(self.bin, JAVA25)
        self.hive_build("build_270/diann-2.7.0")      # nothing to download for DIA-NN

    def tearDown(self):
        for root, dirs, _ in os.walk(self.d):
            for x in dirs:
                os.chmod(os.path.join(root, x), 0o755)
        super().tearDown()

    def run_acquire(self, pin="24.0", extra=None):
        env = {CORE_ENV: os.path.join(self.core, "fragpipe-24.0")}
        env.update(extra or {})
        return self.acquire("hpc", "fragpipe" if pin else "", pin or "", extra_env=env)

    def test_on_hive_the_core_install_is_what_runs_and_is_licensed(self):
        inst = fake_install(self.core)
        t = self.run_acquire()
        self.assertEqual(t["fragpipe"], os.path.join(inst, "bin", "fragpipe"))
        self.assertEqual(t["versions"]["fragpipe"], "24.0")
        ft = t["fragpipe_tools"]
        self.assertEqual((ft["licensed"], ft["ready"], ft["source"], ft["folder"]),
                         (True, True, "core_install", os.path.join(inst, "tools")))
        self.assertEqual(ft["philosopher"]["version"], "5.1.3-RC9")
        self.assertEqual(ft["diann"]["version"], "1.8.2 beta 8")
        self.assertEqual(ft["java"]["major"], 25)
        self.assertFalse(os.path.exists(os.path.join(self.root, "fragpipe")), "nothing downloaded")
        self.assertTrue(any(n.startswith("FragPipe 24.0: ready (tools folder " + ft["folder"])
                            for n in t["notes"]), t["notes"])

    @unittest.skipIf(hasattr(os, "geteuid") and os.geteuid() == 0, "root reads a mode-000 folder")
    def test_an_unreadable_core_install_is_skipped_and_the_licence_note_is_unchanged(self):
        inst = fake_install(self.core)
        os.chmod(inst, 0)
        # the home copy acquire_tools.sh downloaded earlier: no jars
        home = os.path.join(self.root, "fragpipe", "24.0")
        fake_install(home, jars=False)
        t = self.run_acquire()
        self.assertTrue(t["fragpipe"].startswith(home), t["fragpipe"])
        ft = t["fragpipe_tools"]
        self.assertIs(ft["licensed"], False)
        self.assertTrue(any("not readable" in s["why"] for s in ft["skipped"]), ft["skipped"])
        self.assertIn(LICENCE_NOTE.format(ver="24.0"), t["notes"])

    def test_without_python3_nothing_is_claimed(self):
        os.remove(os.path.join(self.bin, "python3"))
        fake_install(self.core)
        t = self.run_acquire()
        self.assertIsNone(t["fragpipe_tools"]["licensed"])
        self.assertIn("python3", t["fragpipe_tools"]["why"])


class DiaTracerJar(unittest.TestCase):
    def test_find_jar_finds_the_capital_T_copy_in_a_subfolder(self):
        import diatracer_parallel as dp
        with tempfile.TemporaryDirectory() as d:
            mine = os.path.join(d, "tools_licensed")
            jar = touch(os.path.join(mine, "diatracer-2.2.1", "diaTracer-2.2.1.jar"))
            old = {k: os.environ.get(k) for k in (CORE_ENV, TOOLS_ENV)}
            os.environ[CORE_ENV] = os.path.join(d, "absent")
            os.environ[TOOLS_ENV] = mine
            try:
                self.assertEqual(dp.find_jar(None), jar)
                tools = os.path.join(d, "tools.json")
                with open(tools, "w") as fh:
                    json.dump({"fragpipe_tools": {"jars": {"diatracer": {"path": jar}}}}, fh)
                os.environ.pop(TOOLS_ENV)
                self.assertEqual(dp.find_jar(None, tools), jar)
            finally:
                for k, v in old.items():
                    if v is None:
                        os.environ.pop(k, None)
                    else:
                        os.environ[k] = v

    def test_the_core_path_is_written_once(self):
        """Rule 3: CORE_FRAGPIPE_INSTALL is the only place the Core install's path is spelled."""
        hits = []
        for name in sorted(os.listdir(SCRIPTS)):
            p = os.path.join(SCRIPTS, name)
            if os.path.isfile(p) and name.endswith((".py", ".sh", ".R")):
                with open(p, errors="replace") as fh:
                    if "proteomics-grp/fragpipe24/fragpipe-24.0" in fh.read():
                        hits.append(name)
        self.assertEqual(hits, ["engine_choice.py"])


if __name__ == "__main__":
    unittest.main()
