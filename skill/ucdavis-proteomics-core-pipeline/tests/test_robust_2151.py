#!/usr/bin/env python3
"""
2.15.1 robustness fixes (audit 2026-11-17, lane 5): everything here is local -- no Docker, no
network, no HIVE. curl, ssh, rsync and docker are stand-ins.

  #2   DIA-NN in Docker (macOS): every absolute path is mounted at the SAME path, --platform
       linux/amd64 is added, and detect_env does not call a stopped Docker "available"
  #3   downloads are checked against the vendor's sha256 (verify_download.sh)
  #4   the Sage on PATH does not override a pin it cannot be shown to match
  #5   HIVE is known by its own marker; another SLURM cluster is not given HIVE's queue
  #7   detect_env.sh on Windows (Git Bash) prints JSON
  #9   --put-skill deletes stale files only inside the skill's own folder on HIVE
  #12  a GitHub failure is said as a network / rate-limit problem, not "no such build"
  #13  a truncated Sage download is not cached
  #18  detect_env JSON is escaped
"""
import hashlib
import io
import json
import os
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
BASH = "/bin/bash"
from job_env import no_cluster_env  # noqa: E402
import test_acquire_tools_versions as tav  # noqa: E402


def _exe(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


class Tmp(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = os.path.realpath(self._t.name)

    def tearDown(self):
        self._t.cleanup()


class DiannDockerCommand(Tmp):
    def dry(self, *args, cwd=None):
        env = dict(os.environ, DIANN_DOCKER_DRYRUN="1")
        r = subprocess.run([BASH, os.path.join(SCRIPTS, "diann_docker.sh"), *args],
                           capture_output=True, text=True, env=env, cwd=cwd or self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        return r.stdout.splitlines()

    def test_cfg_paths_with_spaces_are_mounted_whole_and_never_globbed(self):
        """2.15.1 review minor 13: paths inside a --cfg were split on white space and globbed."""
        lib = os.path.join(self.d, "my libs", "a lib.parquet")
        fasta = os.path.join(self.d, "db*", "x.fasta")
        for p in (lib, fasta):
            os.makedirs(os.path.dirname(p))
            open(p, "w").close()
        os.makedirs(os.path.join(self.d, "dbX"))                # what an unquoted db* would glob to
        cfg = os.path.join(self.d, "p.cfg")
        with open(cfg, "w") as fh:
            fh.write(f"--lib '{lib}'\n--fasta {fasta}\n# --temp /nowhere/else\n--qvalue 0.01\n")
        v = self.dry("img:2.7.0", "--cfg", cfg)
        mounts = [x for x in v if x.startswith("type=bind")]
        self.assertIn(f"type=bind,source={os.path.dirname(lib)},target={os.path.dirname(lib)}",
                      mounts)
        self.assertIn(f"type=bind,source={os.path.dirname(fasta)},target={os.path.dirname(fasta)}",
                      mounts)
        self.assertFalse([m for m in mounts if "dbX" in m or "/my," in m or "nowhere" in m], mounts)

    def test_absolute_paths_are_mounted_at_the_same_path_and_platform_is_set(self):
        raw = os.path.join(self.d, "raw", "a.d")
        fasta = os.path.join(self.d, "db", "x.fasta")
        out = os.path.join(self.d, "out")
        for p in (raw, os.path.dirname(fasta), out):
            os.makedirs(p)
        open(fasta, "w").close()
        v = self.dry("img:2.7.0", "--f", raw, "--fasta", fasta, "--out", os.path.join(out, "report.parquet"),
                     "--temp", os.path.join(out, "quant"), "--threads", "4")
        self.assertEqual(v[:3], ["docker", "run", "--rm"])
        self.assertEqual(v[v.index("--platform") + 1], "linux/amd64")
        mounts = [x for x in v if x.startswith("type=bind")]
        for p in (raw, os.path.dirname(fasta), out):
            self.assertIn(f"type=bind,source={p},target={p}", mounts)
        i = v.index("img:2.7.0")
        self.assertEqual(v[i + 1], "diann-linux")
        self.assertEqual(v[i + 2:i + 5], ["--f", raw, "--fasta"])      # arguments pass through whole
        self.assertEqual(v[v.index("-w") + 1], self.d)

    def test_paths_inside_a_cfg_are_mounted_too_and_the_root_never_is(self):
        lib = os.path.join(self.d, "lib")
        os.makedirs(lib)
        cfg = os.path.join(self.d, "p.cfg")
        with open(cfg, "w") as fh:
            fh.write(f"--lib {lib}/l.speclib\n--threads 4\n")
        v = self.dry("img", "--cfg", cfg, "--out", "/nonexistent_root_xyz/r.parquet")
        mounts = [x for x in v if x.startswith("type=bind")]
        self.assertIn(f"type=bind,source={lib},target={lib}", mounts)
        self.assertNotIn("type=bind,source=/,target=/", mounts)

    def test_acquire_records_the_wrapper_not_a_bare_docker_run(self):
        h = tav.AcquireHarness("setUp")
        h.setUp()
        try:
            tav._exe(os.path.join(h.bin, "docker"), "#!/bin/sh\nexit 0\n")
            t = h.acquire("mac", extra_env={"DIANN_DOCKER_IMAGE": "proteomics-pipeline/diann:2.7.0"})
        finally:
            h.tearDown()
        self.assertIn("diann_docker.sh", t["diann"])
        self.assertIn("proteomics-pipeline/diann:2.7.0", t["diann"])
        self.assertNotIn("$PWD:/data", t["diann"])
        self.assertEqual(t["versions"]["diann"], "2.7.0")


class DetectEnv(Tmp):
    def run_env(self, bin_extra=None, **env):
        b = os.path.join(self.d, "bin")
        os.makedirs(b, exist_ok=True)
        for name, body in (bin_extra or {}).items():
            _exe(os.path.join(b, name), body)
        # ONLY the programs detect_env.sh needs, never /usr/bin whole: a HIVE node has apptainer
        # and java there, and Radiant was then "available" on the simulated Mac (2.15.1 HIVE suite)
        for name in ("uname", "tr", "head", "sed", "dirname", "cat", "perl"):
            if (name not in (bin_extra or {}) and shutil.which(name)
                    and not os.path.lexists(os.path.join(b, name))):
                os.symlink(shutil.which(name), os.path.join(b, name))
        e = no_cluster_env(self.d, base={}, PATH=b, HOME=self.d,
                           HIVE_GROUP_DIR=os.path.join(self.d, "no_group"))
        e.update(env)
        r = subprocess.run([BASH, os.path.join(SCRIPTS, "detect_env.sh")], capture_output=True,
                           text=True, env=e, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        return json.loads(r.stdout)            # must be valid JSON

    def test_windows_git_bash_prints_json_with_no_diann(self):
        j = self.run_env({"uname": "#!/bin/sh\n[ \"$1\" = -m ] && echo x86_64 || echo MINGW64_NT-10.0-22631\n"})
        self.assertFalse(j["engines"]["diann"]["available"])
        self.assertIn("HIVE", j["engines"]["diann"]["note"])

    def test_another_slurm_cluster_is_not_hive_and_needs_an_explicit_queue(self):
        conf = os.path.join(self.d, "cn")
        open(conf, "w").write("farm\n")
        j = self.run_env({"sbatch": "#!/bin/sh\n"}, UCD_HPC_CLUSTER_CONF=conf)
        self.assertEqual((j["platform_class"], j["has_slurm"], j["uc_davis_hive"], j["cluster"]),
                         ("hpc", True, False, "farm"))
        self.assertTrue(j["queue_needs_explicit"])

    def test_hive_is_hive_without_sbatch_on_path(self):
        conf = os.path.join(self.d, "cn")
        open(conf, "w").write('hive\n')
        slurm = os.path.join(self.d, "slurmbin")
        os.makedirs(slurm)
        _exe(os.path.join(slurm, "sbatch"), "#!/bin/sh\n")
        j = self.run_env(UCD_HPC_CLUSTER_CONF=conf, HPCCF_SLURM_BIN=slurm)
        self.assertEqual((j["uc_davis_hive"], j["has_slurm"], j["queue_needs_explicit"]), (True, True, False))

    def test_a_stopped_docker_is_not_available_for_diann_or_radiant(self):
        j = self.run_env({"uname": "#!/bin/sh\n[ \"$1\" = -m ] && echo arm64 || echo Darwin\n",
                          "docker": "#!/bin/sh\nexit 1\n"})
        self.assertTrue(j["has_docker"])
        self.assertFalse(j["docker_running"])
        self.assertFalse(j["engines"]["diann"]["available"])
        self.assertIn("not running", j["engines"]["diann"]["note"])
        self.assertFalse(j["engines"]["radiant"]["available"])

    def test_a_running_docker_is_available_for_diann(self):
        j = self.run_env({"uname": "#!/bin/sh\n[ \"$1\" = -m ] && echo arm64 || echo Darwin\n",
                          "docker": "#!/bin/sh\nexit 0\n"})
        self.assertTrue(j["engines"]["diann"]["available"])

    def test_text_with_quotes_and_backslashes_stays_valid_json(self):
        # a Java banner with a backslash: the old interpolated JSON broke on it
        j = self.run_env({"java": "#!/bin/sh\necho 'openjdk version \\1.8 C:\\jdk'\n"})
        self.assertIn("jdk", j["has_java"])


class SlurmQueueNotHive(Tmp):
    def queue(self, conf_text, *args, assoc="a1|p1|q1\na2|p2|q2\n", group_dir=None):
        b = os.path.join(self.d, "bin")
        os.makedirs(b, exist_ok=True)
        _exe(os.path.join(b, "sbatch"), "#!/bin/sh\n")
        _exe(os.path.join(b, "sacctmgr"), f"#!/bin/sh\ncat <<'E'\n{assoc}E\n")
        conf = os.path.join(self.d, "cn")
        if conf_text is not None:
            open(conf, "w").write(conf_text)
        e = no_cluster_env(self.d, base={}, PATH=f"{b}:/usr/bin:/bin", USER="tester",
                           UCD_HPC_CLUSTER_CONF=conf,
                           HIVE_GROUP_DIR=group_dir or os.path.join(self.d, "no_group"))
        code = ("import sys; sys.path.insert(0, %r); import run_search as r; "
                "print(r.slurm_queue(*%r))" % (SCRIPTS, args))
        return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                              env=e, cwd=self.d, timeout=60)

    def test_cluster_is_hive_by_its_own_marker(self):
        code = ("import sys; sys.path.insert(0, %r); import run_search as r, os; "
                "print(r.cluster_is_hive())" % SCRIPTS)
        for text, want in (("hive\n", "True"), ("farm\n", "False"), ("Franklin", "False"),
                           ('export CLUSTER_NAME="hive"\n', "True")):
            conf = os.path.join(self.d, "cn")
            open(conf, "w").write(text)
            e = dict(os.environ, UCD_HPC_CLUSTER_CONF=conf, HIVE_GROUP_DIR="/nonexistent")
            r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=e)
            self.assertEqual(r.stdout.strip(), want, text)

    def test_a_non_hive_cluster_with_several_associations_is_not_given_hives_queue(self):
        r = self.queue("farm\n")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not HIVE", r.stderr)
        self.assertIn("--partition", r.stderr)
        self.assertNotIn("publicgrp", r.stdout)

    def test_a_non_hive_cluster_with_unreadable_associations_asks_too(self):
        r = self.queue("farm\n", assoc="")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not HIVE", r.stderr)

    def test_a_single_association_is_used_as_it_is(self):
        r = self.queue("farm\n", assoc="onlyacct|onlypart|onlyqos\n")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stdout.strip(), "('onlypart', 'onlyacct', 'onlyqos')")

    def test_an_explicit_pair_is_honoured_on_any_cluster(self):
        r = self.queue("farm\n", "p1", "a1")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("'p1'", r.stdout)

    def test_on_hive_the_old_rule_stands(self):
        r = self.queue("hive\n", assoc="publicgrp|low|publicgrp-low-qos\n")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("publicgrp", r.stdout)


class VerifyDownload(Tmp):
    def verify(self, table_lines, key, content=b"hello"):
        tbl = os.path.join(self.d, "t.tsv")
        with open(tbl, "w") as fh:
            fh.write("\n".join(table_lines) + "\n")
        f = os.path.join(self.d, "f.bin")
        open(f, "wb").write(content)
        r = subprocess.run([BASH, "-c", f'. "{SCRIPTS}/verify_download.sh"; verify_download {key} "{f}" thing; echo "rc=$?"'],
                           capture_output=True, text=True, env=dict(os.environ, PINNED_CHECKSUMS=tbl))
        return r.stdout.strip(), r.stderr

    def test_match_mismatch_and_not_pinned(self):
        good = hashlib.sha256(b"hello").hexdigest()
        out, err = self.verify([f"k1 {good} source"], "k1")
        self.assertEqual(out, "rc=0")
        out, err = self.verify([f"k1 {'0' * 64} source"], "k1")
        self.assertEqual(out, "rc=1")
        self.assertIn("NOT used", err)
        out, err = self.verify([f"k1 {good} source"], "other")
        self.assertEqual(out, "rc=2")
        self.assertIn("NOT verified", err)

    def test_the_shipped_table_has_a_value_for_every_pinned_engine_that_publishes_one(self):
        with open(os.path.join(SCRIPTS, "pinned_checksums.txt")) as fh:
            keys = {ln.split()[0]: ln.split()[1] for ln in fh if ln.strip() and not ln.startswith("#")}
        for k in ("fragpipe-24.0-linux", "diann-2.7.0-linux", "dotnet-install-5147e32",
                  "micromamba-2.10.0-0-osx-arm64", "micromamba-2.10.0-0-linux-64"):
            self.assertIn(k, keys)
        for v in keys.values():
            self.assertRegex(v, r"^[0-9a-f]{64}$")
        self.assertNotIn("sage-0.14.7-x86_64-unknown-linux-gnu", keys)   # Sage publishes none


class AcquireDownloads(tav.AcquireHarness):
    def table(self, line):
        p = os.path.join(self.d, "pins.tsv")
        open(p, "w").write(line + "\n")
        return p

    def fp_fixtures(self):
        with open(os.path.join(self.fix, "fragpipe_releases.json"), "w") as fh:
            fh.write('[{"assets":[{"browser_download_url":"https://github.com/Nesvilab/FragPipe/releases/download/24.0/FragPipe-24.0-linux.zip"}]}]')
            fh.write("\n")
        import zipfile
        with zipfile.ZipFile(os.path.join(self.fix, "fragpipe.zip"), "w") as z:
            info = zipfile.ZipInfo("fragpipe-24.0/bin/fragpipe")
            info.external_attr = 0o755 << 16
            z.writestr(info, "#!/bin/sh\n")

    def test_a_fragpipe_zip_that_does_not_match_its_published_sha_is_not_installed(self):
        self.fp_fixtures()
        # the releases json for the tag endpoint: the stub serves one file for any FragPipe API call
        t = self.acquire("linux", "fragpipe", "24.0", extra_env={
            "PINNED_CHECKSUMS": self.table(f"fragpipe-24.0-linux {'0' * 64} test")})
        self.assertIsNone(t["fragpipe"])
        self.assertTrue(any("did NOT match" in n for n in t["notes"]), t["notes"])
        self.assertFalse(os.path.exists(os.path.join(self.root, "fragpipe", "24.0", "fragpipe.zip")))

    def test_a_fragpipe_zip_that_matches_is_installed(self):
        self.fp_fixtures()
        with open(os.path.join(self.fix, "fragpipe.zip"), "rb") as fh:
            sha = hashlib.sha256(fh.read()).hexdigest()
        t = self.acquire("linux", "fragpipe", "24.0", extra_env={
            "PINNED_CHECKSUMS": self.table(f"fragpipe-24.0-linux {sha} test")})
        self.assertTrue(t["fragpipe"] and t["fragpipe"].endswith("bin/fragpipe"), t["fragpipe"])

    def test_diann_zip_mismatch_is_not_installed(self):
        t = self.acquire("linux", "diann", "2.7.0", extra_env={
            "PINNED_CHECKSUMS": self.table(f"diann-2.7.0-linux {'1' * 64} test")})
        self.assertIsNone(t["diann"])
        self.assertTrue(any("did NOT match" in n for n in t["notes"]), t["notes"])

    def test_unpinned_download_is_said_to_be_not_verified(self):
        t = self.acquire("linux", "diann", "2.7.0")
        self.assertTrue(t["diann"])
        self.assertTrue(any("NOT verified" in n for n in t["notes"]), t["notes"])


class SageAcquire(tav.AcquireHarness):
    def sage_tar(self, truncated=False):
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tf:
            top = "sage-v0.14.7-x86_64-unknown-linux-gnu"
            ti = tarfile.TarInfo(top)
            ti.type = tarfile.DIRTYPE
            tf.addfile(ti)
            body = b"#!/bin/sh\necho sage\n" + b"x" * 20000
            ti = tarfile.TarInfo(f"{top}/sage")
            ti.size = len(body)
            tf.addfile(ti, io.BytesIO(body))
        data = buf.getvalue()
        return data[: len(data) // 2] if truncated else data

    def stub(self, tar_bytes):
        with open(os.path.join(self.fix, "sage.tar.gz"), "wb") as fh:
            fh.write(tar_bytes)
        with open(os.path.join(self.fix, "sage_release.json"), "w") as fh:
            fh.write('{"assets":[{"browser_download_url":"https://github.com/lazear/sage/releases/download/v0.14.7/sage-v0.14.7-x86_64-unknown-linux-gnu.tar.gz"}]}\n')
        stub = open(os.path.join(self.bin, "curl")).read()
        stub = stub.replace("  *) exit 22 ;;\nesac", '''  *api.github.com/repos/lazear/sage/*) emit < "$FIXTURES/sage_release.json" ;;
  *sage-v*.tar.gz) cp "$FIXTURES/sage.tar.gz" "$out" ;;
  *) exit 22 ;;
esac''')
        tav._exe(os.path.join(self.bin, "curl"), stub)

    def _on_linux(self):
        # the harness links the real uname; force linux x86_64
        tav._exe(os.path.join(self.bin, "uname"), "#!/bin/sh\n[ \"$1\" = -m ] && echo x86_64 || echo Linux\n")

    def test_a_truncated_download_is_not_cached_as_sage(self):
        self._on_linux()
        self.stub(self.sage_tar(truncated=True))
        t = self.acquire("linux", "sage", "0.14.7")
        self.assertIsNone(t["sage"])
        d = os.path.join(self.root, "sage", "0.14.7")
        self.assertFalse(os.path.exists(os.path.join(d, "sage")))
        self.assertTrue(any("damaged" in n or "cut short" in n for n in t["notes"]), t["notes"])

    def test_a_whole_download_is_installed_and_says_not_verified(self):
        self._on_linux()
        self.stub(self.sage_tar())
        t = self.acquire("linux", "sage", "0.14.7")
        self.assertTrue(t["sage"].endswith("sage/0.14.7/sage"), t["sage"])
        self.assertEqual(t["versions"]["sage"], "0.14.7")
        self.assertTrue(any("NOT verified" in n for n in t["notes"]), t["notes"])

    def test_a_sage_on_path_that_is_not_the_pin_does_not_override_it(self):
        self._on_linux()
        self.stub(self.sage_tar())
        tav._exe(os.path.join(self.bin, "sage"), "#!/bin/sh\necho 'sage 0.13.0'\n")
        t = self.acquire("linux", "sage", "0.14.7")
        self.assertTrue(t["sage"].endswith("sage/0.14.7/sage"), t["sage"])
        self.assertTrue(any("NOT used" in n and "0.13.0" in n for n in t["notes"]), t["notes"])

    def test_a_sage_on_path_that_is_the_pin_is_used_and_records_the_version(self):
        self._on_linux()
        tav._exe(os.path.join(self.bin, "sage"), "#!/bin/sh\necho 'sage 0.14.7'\n")
        t = self.acquire("linux", "sage", "0.14.7")
        self.assertEqual(t["sage"], os.path.join(self.bin, "sage"))
        self.assertEqual(t["versions"]["sage"], "0.14.7")

    def test_an_unpinned_run_still_takes_the_sage_on_path(self):
        tav._exe(os.path.join(self.bin, "sage"), "#!/bin/sh\necho 'sage 0.13.0'\n")
        t = self.acquire("linux")
        self.assertEqual(t["versions"]["sage"], "env")


class GitHubFailureWording(tav.AcquireHarness):
    def test_an_unreachable_github_is_not_reported_as_no_such_build(self):
        tav._exe(os.path.join(self.bin, "curl"), "#!/bin/sh\nexit 6\n")
        t = self.acquire("linux", "diann", "2.7.0")
        text = " ".join(t["notes"])
        self.assertIn("could not be reached", text)
        self.assertNotIn("no Academia Linux build named", text)

    def test_a_rate_limit_is_said_as_one(self):
        tav._exe(os.path.join(self.bin, "curl"),
                 '#!/bin/sh\nout=""; while [ $# -gt 0 ]; do [ "$1" = -o ] && out="$2"; shift; done\n'
                 '[ -n "$out" ] && echo \'{"message":"rate limit"}\' > "$out"; printf 403\n')
        t = self.acquire("linux", "diann", "2.7.0")
        text = " ".join(t["notes"])
        self.assertIn("rate limit", text)
        self.assertNotIn("no Academia Linux build named", text)

    def test_a_missing_version_is_still_said_as_missing(self):
        t = self.acquire("linux", "diann", "9.9.9")
        self.assertTrue(any("no Academia Linux build named '9.9.9'" in n for n in t["notes"]), t["notes"])


class PutSkill(Tmp):
    """--put-skill with a fake ssh (runs the remote command here) and a logging rsync."""

    def run_put(self, with_rsync=True):
        b = os.path.join(self.d, "bin")
        os.makedirs(b)
        for t in ("tar", "mkdir", "cp", "rm", "mktemp", "dirname", "basename", "cat", "ls", "tr", "uname",
                  "mv", "rmdir", "sed", "grep", "head", "awk", "sort", "wc", "chmod", "date", "git",
                  "printf", "id"):
            real = shutil.which(t)
            if real:
                os.symlink(real, os.path.join(b, t))
        os.symlink(BASH, os.path.join(b, "bash"))
        logs = os.path.join(self.d, "logs")
        os.makedirs(logs)
        _exe(os.path.join(b, "ssh"), '#!/bin/bash\nfor last; do :; done\nexec bash -c "$last"\n')
        if with_rsync:
            _exe(os.path.join(b, "rsync"), '#!/bin/bash\nfor a in "$@"; do printf "%s\\n" "$a"; done >> "$LOGS/rsync.log"\necho --end-- >> "$LOGS/rsync.log"\n')
        root = os.path.join(self.d, "skill")
        os.makedirs(os.path.join(root, ".claude-plugin"))
        os.makedirs(os.path.join(root, "scripts"))
        open(os.path.join(root, ".claude-plugin", "plugin.json"), "w").write('{\n  "name": "ucdavis-proteomics-core-pipeline",\n  "version": "2.15.1"\n}\n')
        shutil.copy(os.path.join(SCRIPTS, "skill_version.sh"), os.path.join(root, "scripts", "skill_version.sh"))
        shutil.copy(os.path.join(SCRIPTS, "hive_exec.sh"), os.path.join(root, "scripts", "hive_exec.sh"))
        key = os.path.join(self.d, "k")
        open(key, "w").close()
        home = os.path.join(self.d, "home")
        os.makedirs(os.path.join(home, "proteomics-pipeline", ".claude-plugin"))
        open(os.path.join(home, "proteomics-pipeline", ".claude-plugin", "plugin.json"), "w").write('{\n  "version": "2.15.0"\n}\n')
        env = {"PATH": b, "HOME": home, "LOGS": logs, "HIVE_USER": "tester", "HIVE_KEY": key,
               "HIVE_ENV_FILE": os.path.join(self.d, "no.env"), "HIVE_SSH_MUX": "0",
               "HIVE_PUT_FORCE": "1", "TMPDIR": self.d}
        script = os.path.join(root, "scripts", "hive_exec.sh")
        r = subprocess.run([BASH, script, "--put-skill"], capture_output=True, text=True, env=env, timeout=60)
        return r, logs, home

    def test_rsync_deletes_only_inside_the_skills_own_folders_and_records_the_put(self):
        r, logs, home = self.run_put()
        self.assertEqual(r.returncode, 0, r.stderr)
        log = open(os.path.join(logs, "rsync.log")).read().splitlines()
        calls, cur = [], []
        for ln in log:
            if ln == "--end--":
                calls.append(cur); cur = []
            else:
                cur.append(ln)
        self.assertEqual(len(calls), 2)
        dests = []
        for c in calls:
            self.assertIn("--delete", c)
            dests.append(c[-1])
        self.assertEqual(dests, ["tester@hive.hpc.ucdavis.edu:~/proteomics-pipeline/scripts/",
                                 "tester@hive.hpc.ucdavis.edu:~/proteomics-pipeline/.claude-plugin/"])
        for c in calls:
            self.assertTrue(c[-2].endswith("/"), c)          # the CONTENTS of the skill's own folder
        rec = open(os.path.join(home, "proteomics-pipeline", ".put_record")).read()
        self.assertIn("2.15.0 -> 2.15.1", rec)
        self.assertIn("tester", rec)


class DetectAcquisitionPaths(Tmp):
    def run_da(self, *args):
        return subprocess.run([sys.executable, "-I", os.path.join(SCRIPTS, "detect_acquisition.py"), *args],
                              capture_output=True, text=True, cwd=self.d, timeout=60)

    def test_a_missing_path_is_an_error_not_an_unknown_run(self):
        r = self.run_da(os.path.join(self.d, "nope.raw"))
        self.assertEqual(r.returncode, 2)
        self.assertIn("not found", r.stderr)
        self.assertEqual(r.stdout, "")

    def test_an_empty_folder_is_an_error(self):
        os.makedirs(os.path.join(self.d, "empty"))
        r = self.run_da(os.path.join(self.d, "empty"))
        self.assertEqual(r.returncode, 2)
        self.assertIn("no .d, .raw, .mzML or .wiff runs", r.stderr)

    def test_a_folder_is_read_as_the_runs_inside_it(self):
        f = os.path.join(self.d, "runs")
        os.makedirs(os.path.join(f, "a.d"))
        open(os.path.join(f, "notes.txt"), "w").close()
        r = self.run_da(f)
        self.assertEqual(r.returncode, 0, r.stderr)
        files = [x["file"] for x in json.loads(r.stdout)["files"]]
        self.assertEqual([os.path.basename(x) for x in files], ["a.d"])

    def test_help_exits_zero(self):
        self.assertEqual(self.run_da("--help").returncode, 0)


class StopHookIsHonest(unittest.TestCase):
    def test_hooks_json_says_a_hive_session_is_not_checked(self):
        with open(os.path.join(os.path.dirname(HERE), "hooks", "hooks.json")) as fh:
            d = json.load(fh)["description"]
        self.assertIn("NOT check a session whose folder is on HIVE", d)
        self.assertIn("session.py status", d)


class EnvironmentIsPinned(unittest.TestCase):
    """setup.sh (finding 4 and 3): versions are pinned and micromamba is not `latest`."""

    def setUp(self):
        with open(os.path.join(SCRIPTS, "setup.sh")) as fh:
            self.src = fh.read()

    def test_every_analysis_package_named_in_the_finding_carries_a_pin(self):
        for pkg in ("bioconductor-limma=", "r-statmod=", "r-arrow=", "pyarrow=", "sage-proteomics="):
            self.assertIn(f'"{pkg}$PIN_', self.src, pkg)
        self.assertRegex(self.src, r'PIN_LIMMA="\$\{PIN_LIMMA:-3\.68\.5\}"')
        self.assertIn('LIMPA_PIN="$PIN_LIMPA"', self.src)

    def test_micromamba_is_a_pinned_release_checked_against_its_sha(self):
        code = "\n".join(ln for ln in self.src.splitlines() if not ln.lstrip().startswith("#"))
        self.assertNotIn("curl -Ls", code)
        self.assertNotIn("micro.mamba.pm/api", code)
        self.assertIn("verify_download", self.src)
        self.assertRegex(self.src, r'PIN_MICROMAMBA="\$\{PIN_MICROMAMBA:-\d+\.\d+\.\d+-\d+\}"')

    def test_the_pinned_sage_release_matches_resolve_defaults(self):
        sys.path.insert(0, SCRIPTS)
        import resolve_defaults as rd
        pins = getattr(rd, "ENGINE_PINS", None) or getattr(rd, "PINS", None)
        if isinstance(pins, dict) and "sage" in pins:
            self.assertEqual(pins["sage"], "0.14.7")
        self.assertIn("0.14.7", self.src)


if __name__ == "__main__":
    unittest.main()
