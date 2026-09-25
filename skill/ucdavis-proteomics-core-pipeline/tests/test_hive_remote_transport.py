#!/usr/bin/env python3
"""
hive_remote from a Windows laptop: what one real session (gabrig, 2026-09-23, Windows 11 +
Git Bash) ran into, and what now stops it.

  1. The raw folder was T:\\Data\\lab\\service\\..., and T: is \\\\128.120.208.24\\proteomics --
     the Flinders share HIVE mounts at /nfs/lssc0/flinders/proteomics. Nothing mapped the one
     to the other; 7.1 GB that was already on HIVE was uploaded (~20 min). hive_path.sh maps
     a local path to its HIVE path and verifies it, and `hive_exec.sh --put` refuses to
     upload a source it verifies.
  2. Git Bash has no rsync, so --put/--get died. They fall back to scp -- keeping rsync's
     meaning, which scp -r does not have on its own.
  3. First contact failed "Host key verification failed"; the default ssh user was the AD
     login `AD3+gabrig` ("Permission denied"); python3 was the Microsoft Store stub.
     check_access.sh now says which of these it is, instead of hive_ssh: "failed".

No test contacts a real host: ssh, scp, rsync, ssh-keygen and ssh-keyscan are recorders on a
PATH that holds nothing else but a few linked coreutils, `uname`/`net`/`mount`/`cygpath` are
fakes, and HIVE itself is a directory tree reached through HIVE_EXEC.
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
HIVE_PATH = os.path.join(SCRIPTS, "hive_path.sh")
HIVE_EXEC = os.path.join(SCRIPTS, "hive_exec.sh")
CHECK_ACCESS = os.path.join(SCRIPTS, "check_access.sh")

COREUTILS = ("bash", "tr", "awk", "dirname", "basename", "ls", "sort", "wc", "stat", "grep",
             "head", "tail", "cat", "mkdir", "mv", "rmdir", "rm", "mktemp", "sed")
HOST = "hive.hpc.ucdavis.edu"
FLINDERS = "/nfs/lssc0/flinders/proteomics"
PUBLISHED_ED25519 = "SHA256:b5nv86Ciaqg1yrUVai6bZ0Hk4IpzAFLWtIPDBdacbQM"

NET_USE = (
    "New connections will be remembered.\r\n\r\n\r\n"
    "Status       Local     Remote                    Network\r\n\r\n"
    "-------------------------------------------------------------------------------\r\n"
    "OK           T:        \\\\128.120.208.24\\proteomics\r\n"
    "                                                Microsoft Windows Network\r\n"
    "Unavailable  R:        \\\\files.example.edu\\bioinfo  Microsoft Windows Network\r\n"
    "The command completed successfully.\r\n\r\n"
)

# Every fake is a small bash script. RECORD appends one call's argv, one arg per line.
RECORD = 'for a in "$@"; do printf "%s\\n" "$a"; done >> "$LOGS/{name}.log"; echo "--end--" >> "$LOGS/{name}.log"\n'
SHIMS = {
    "uname": 'echo "${FAKE_UNAME:-Darwin}"\n',
    "mount": 'cat "$FAKE_MOUNT" 2>/dev/null\n',
    "net": '[ -f "$FAKE_NET" ] || exit 1\ncat "$FAKE_NET"\n',
    "powershell.exe": RECORD.format(name="powershell") + 'printf "%s\\r\\n" "$FAKE_PS_ROOT"\n',
    # cygpath -u <Windows path>: the share's local copy lives under $SHAREROOT
    "cygpath": ('[ -n "${FAKE_CYGPATH_OUT:-}" ] && { echo "$FAKE_CYGPATH_OUT"; exit; }\n'
                'p="${2//\\\\//}"\n'
                'case "$p" in [A-Za-z]:*) p="${p#?:}" ;; //*) p="${p#//}"; p="${p#*/}"; p="/${p#*/}" ;; esac\n'
                'echo "$SHAREROOT$p"\n'),
    "ssh": RECORD.format(name="ssh") + (
        'case "${FAKE_SSH:-ok}" in\n'
        '  ok) echo HAS_SBATCH; echo HAS_GRP ;;\n'
        '  old_ls_d) echo HAS_SBATCH; echo /quobyte/proteomics-grp ;;\n'
        '  hostkey) echo "Warning: Permanently added x" >&2; echo "Host key verification failed." >&2; exit 255 ;;\n'
        '  denied) echo "UC Davis HPC -- authorised use only" >&2\n'
        '          echo "tester@hive.hpc.ucdavis.edu: Permission denied (publickey)." >&2; exit 255 ;;\n'
        '  throttled) echo "kex_exchange_identification: read: Operation timed out" >&2; exit 255 ;;\n'
        '  other) echo "UC Davis HPC -- authorised use only" >&2\n'
        '         echo "Received disconnect from 1.2.3.4 port 22:2: Too many authentication failures" >&2\n'
        # OpenSSH always follows a Received disconnect with this reasonless line (packet.c)
        '         echo "Disconnected from 1.2.3.4 port 22" >&2; exit 255 ;;\n'
        'esac\n'),
    "scp": RECORD.format(name="scp") + (
        # --get into a directory: leave a file named after the remote source, as scp would
        'if [ -n "${FAKE_SCP_CREATE:-}" ]; then\n'
        '  for last; do :; done\n'
        '  [ -d "$last" ] && echo data > "$last/${FAKE_SCP_CREATE}"\n'
        'fi\n'),
    "rsync": RECORD.format(name="rsync"),
    "ssh-keygen": RECORD.format(name="ssh-keygen") + (
        'case "$*" in\n'
        '  *-lf*) cat >/dev/null; echo "256 ${FAKE_SCAN_FP} hive.hpc.ucdavis.edu (ED25519)" ;;\n'
        '  *-F*/etc/ssh*) exit 1 ;;\n'
        '  *-F*) [ -n "${FAKE_KNOWN_FP:-}" ] || exit 1\n'
        '        echo "# Host hive.hpc.ucdavis.edu found: line 1"\n'
        '        echo "hive.hpc.ucdavis.edu ED25519 ${FAKE_KNOWN_FP}" ;;\n'
        'esac\n'),
    "ssh-keyscan": RECORD.format(name="ssh-keyscan") + 'echo "hive.hpc.ucdavis.edu ssh-ed25519 AAAAC3fake"\n',
}

# HIVE_EXEC stand-in: runs the remote script against $HIVEROOT instead of /, and maps the
# paths it prints back, so the script under test sees exactly what HIVE would answer.
FAKE_HIVE = (
    'echo call >> "$LOGS/hive.log"\n'
    's="${1//\\/nfs\\/lssc0\\//$HIVEROOT/nfs/lssc0/}"\n'
    's="${s//\\/quobyte\\//$HIVEROOT/quobyte/}"\n'
    'bash -c "$s" 2>&1 | while IFS= read -r l; do printf "%s\\n" "${l//$HIVEROOT/}"; done\n'
)


def _write_exe(path, body):
    with open(path, "w") as fh:
        fh.write("#!/bin/bash\n" + body)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


class Harness(unittest.TestCase):
    SKIP_SHIMS = ()

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        # realpath: macOS's /var is /private/var, and `cd -P` in the scripts resolves it
        self.d = os.path.realpath(self._tmp.name)
        self.bin = self._mk("bin")
        self.logs = self._mk("logs")
        self.home = self._mk("home")
        for tool in COREUTILS:
            real = shutil.which(tool)
            if not real:
                self.skipTest(f"{tool} is not installed")
            os.symlink(real, os.path.join(self.bin, tool))
        for name, body in SHIMS.items():
            if name not in self.SKIP_SHIMS:
                _write_exe(os.path.join(self.bin, name), body)
        self.fake_hive = os.path.join(self.d, "fake_hive.sh")
        _write_exe(self.fake_hive, FAKE_HIVE)
        self.key = os.path.join(self.d, "id_test")
        open(self.key, "w").close()
        # what the laptop sees on the share, and what HIVE has at /nfs/lssc0/flinders/proteomics
        self.shareroot = self._mk("share")
        self.hiveroot = self._mk("hive")
        self.service = "Data/lab/service"
        for root in (self.shareroot, os.path.join(self.hiveroot, FLINDERS.lstrip("/"))):
            run = os.path.join(root, self.service, "P1")
            os.makedirs(os.path.join(run, "c.d"))
            for f in ("a.raw", "b.raw"):
                with open(os.path.join(run, f), "w") as fh:
                    fh.write("x" * 10)
            with open(os.path.join(root, self.service, "one.raw"), "w") as fh:
                fh.write("y" * 1234)
        # an instrument still writing: HIVE may hold MORE than the laptop listed
        open(os.path.join(self.hiveroot, FLINDERS.lstrip("/"), self.service, "P1", "new.raw"), "w").close()
        self.mount_file = os.path.join(self.d, "mount.txt")
        self.set_mount("/dev/disk3s1s1 on / (apfs, sealed, local, read-only, journaled)\n")
        self.net_file = os.path.join(self.d, "net.txt")
        with open(self.net_file, "w", newline="") as fh:
            fh.write(NET_USE)

    def tearDown(self):
        self._tmp.cleanup()

    def _mk(self, *parts):
        p = os.path.join(self.d, *parts)
        os.makedirs(p, exist_ok=True)
        return p

    def set_mount(self, text):
        with open(self.mount_file, "w") as fh:
            fh.write(text)

    def env(self, **extra):
        e = {"PATH": self.bin, "HOME": self.home, "LOGS": self.logs, "HIVEROOT": self.hiveroot,
             "SHAREROOT": self.shareroot, "FAKE_MOUNT": self.mount_file, "FAKE_NET": self.net_file,
             "HIVE_EXEC": f"bash {self.fake_hive}", "HIVE_USER": "tester", "HIVE_KEY": self.key,
             "HIVE_ENV_FILE": os.path.join(self.d, "no.env"), "TMPDIR": self._mk("tmp")}
        e.update(extra)
        return e

    def run_script(self, script, *args, **env):
        # job_env: not a search job (copies of the skill's own hive_exec.sh / hive_path.sh)
        return subprocess.run(["bash", script, *args], capture_output=True, text=True,
                              env=self.env(**env), cwd=self.d, timeout=60)

    def _log(self, name):
        p = os.path.join(self.logs, f"{name}.log")
        if not os.path.exists(p):
            return ""
        with open(p) as fh:
            return fh.read()

    def calls(self, name):
        """argv of every call a recorder saw, in order"""
        out, cur = [], []
        for line in self._log(name).split("\n"):
            if line == "--end--":
                out.append(cur)
                cur = []
            elif line or cur:
                cur.append(line)
        return out

    def hive_calls(self):
        return len(self._log("hive").split())


class HivePathWindowsTests(Harness):
    WIN = {"FAKE_UNAME": "MINGW64_NT-10.0-22631"}

    def resolve(self, path, **env):
        r = self.run_script(HIVE_PATH, path, **{**self.WIN, **env})
        return r, json.loads(r.stdout)

    def test_a_mapped_drive_resolves_to_the_flinders_mount_and_is_verified(self):
        """gabrig's exact path shape. One ssh call, and HIVE having an extra file the laptop
        has not seen yet (an instrument still writing) does not break the match."""
        for form in ("T:\\Data\\lab\\service\\P1", "T:/Data/lab/service/P1"):
            with self.subTest(form=form):
                r, j = self.resolve(form)
                self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
                self.assertEqual(j["unc_or_mount_source"], "\\\\128.120.208.24\\proteomics")
                self.assertEqual((j["server"], j["share"], j["rest"]),
                                 ("128.120.208.24", "proteomics", "Data/lab/service/P1"))
                self.assertEqual(j["candidates"], [f"{FLINDERS}/Data/lab/service/P1"])
                self.assertIs(j["verified"], True)
                self.assertEqual(j["hive_path"], f"{FLINDERS}/Data/lab/service/P1")
        self.assertEqual(self.hive_calls(), 2)          # one per resolve, never one per candidate

    def test_unc_paths_in_both_spellings(self):
        for form in ("\\\\128.120.208.24\\proteomics\\Data\\lab\\service\\P1",
                     "//128.120.208.24/proteomics/Data/lab/service/P1"):
            with self.subTest(form=form):
                r, j = self.resolve(form)
                self.assertEqual((j["server"], j["share"], j["rest"]),
                                 ("128.120.208.24", "proteomics", "Data/lab/service/P1"))
                self.assertIs(j["verified"], True, j["how"])

    def test_git_bash_drive_spelling_is_parsed_as_the_drive(self):
        """/t/Data/... is how Git Bash shows T:\\Data\\... . It is not readable on this test
        machine, so it cannot be verified -- and must not claim to be."""
        r, j = self.resolve("/t/Data/lab/service/P1")
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertEqual((j["share"], j["rest"]), ("proteomics", "Data/lab/service/P1"))
        self.assertEqual(j["candidates"], [f"{FLINDERS}/Data/lab/service/P1"])
        self.assertIs(j["verified"], False)
        self.assertIsNone(j["hive_path"])
        self.assertEqual(self.hive_calls(), 0)

    def test_a_file_is_verified_by_size(self):
        r, j = self.resolve("T:\\Data\\lab\\service\\one.raw")
        self.assertIs(j["verified"], True, j["how"])
        self.assertIn("1234 bytes", j["how"])
        with open(os.path.join(self.hiveroot, FLINDERS.lstrip("/"), self.service, "one.raw"), "a") as fh:
            fh.write("more")
        r, j = self.resolve("T:\\Data\\lab\\service\\one.raw")
        self.assertEqual(r.returncode, 1)
        self.assertIs(j["verified"], False)
        self.assertIsNone(j["hive_path"])

    def test_a_directory_missing_names_on_hive_is_not_verified(self):
        open(os.path.join(self.shareroot, self.service, "P1", "only_here.raw"), "w").close()
        r, j = self.resolve("T:\\Data\\lab\\service\\P1")
        self.assertIs(j["verified"], False)
        self.assertIn("1 of 4 names missing", j["how"])

    def test_case_typed_differently_on_windows_still_finds_the_hive_folder(self):
        """SMB ignores case, so t:\\data\\LAB works on the laptop; HIVE's NFS does not."""
        r, j = self.resolve("t:\\data\\LAB\\service\\P1",
                            FAKE_CYGPATH_OUT=os.path.join(self.shareroot, self.service, "P1"))
        self.assertIs(j["verified"], True, j["how"])
        # A case-insensitive test filesystem (macOS APFS) answers `data` itself, so only a
        # case-sensitive one -- like HIVE's NFS -- shows the fold happened.
        if not os.path.exists(os.path.join(self.hiveroot, FLINDERS.lstrip("/"), "data")):
            self.assertEqual(j["hive_path"], f"{FLINDERS}/Data/lab/service/P1")
        self.assertEqual(j["hive_path"].lower(), f"{FLINDERS}/Data/lab/service/P1".lower())

    def test_a_local_disk_is_not_a_share_and_makes_no_ssh_call(self):
        r, j = self.resolve("C:\\Users\\gabrig\\raw")
        self.assertEqual(r.returncode, 3, r.stdout + r.stderr)
        self.assertEqual(j["candidates"], [])
        self.assertEqual(self.hive_calls(), 0)

    def test_powershell_is_the_fallback_when_there_is_no_net(self):
        os.remove(os.path.join(self.bin, "net"))
        r, j = self.resolve("T:\\Data\\lab\\service\\P1", FAKE_PS_ROOT="\\\\128.120.208.24\\proteomics")
        self.assertEqual(j["share"], "proteomics")
        self.assertIs(j["verified"], True, j["how"])
        self.assertEqual(len(self.calls("powershell")), 1)

    def test_an_unknown_share_tries_flinders_then_quobyte_in_one_call(self):
        bio = os.path.join(self.hiveroot, "quobyte", "bioinfo", "runs")
        os.makedirs(bio)
        local = self._mk("share", "runs")
        for p in (bio, local):
            open(os.path.join(p, "r1.raw"), "w").close()
        r, j = self.resolve("R:\\runs")
        self.assertEqual(j["candidates"], ["/nfs/lssc0/flinders/bioinfo/runs", "/quobyte/bioinfo/runs"])
        self.assertIs(j["verified"], True, j["how"])
        self.assertEqual(j["hive_path"], "/quobyte/bioinfo/runs")
        self.assertEqual(self.hive_calls(), 1)


class HivePathMountTests(Harness):
    def test_a_macos_smb_mount_of_the_flinders_share(self):
        vol = self._mk("Volumes", "proteomics")
        shutil.rmtree(vol)
        shutil.copytree(self.shareroot, vol)
        self.set_mount("/dev/disk3s1s1 on / (apfs, sealed, local, read-only, journaled)\n"
                       f"//gc-prot-core-user@128.120.208.24/proteomics on {vol} (smbfs, nodev, nosuid, mounted by x)\n")
        r = self.run_script(HIVE_PATH, os.path.join(vol, self.service, "P1"))
        j = json.loads(r.stdout)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(j["unc_or_mount_source"], "//gc-prot-core-user@128.120.208.24/proteomics")
        self.assertEqual(j["hive_path"], f"{FLINDERS}/Data/lab/service/P1")

    def test_a_relative_path_inside_the_mount(self):
        vol = self._mk("Volumes", "proteomics")
        shutil.rmtree(vol)
        shutil.copytree(self.shareroot, vol)
        self.set_mount(f"//u@128.120.208.24/proteomics on {vol} (smbfs, nodev)\n")
        r = subprocess.run(["bash", HIVE_PATH, "P1"], capture_output=True, text=True,
                           env=self.env(), cwd=os.path.join(vol, self.service))
        self.assertIs(json.loads(r.stdout)["verified"], True, r.stdout)

    def test_proteomics_grp_maps_to_quobyte(self):
        vol = self._mk("Volumes", "proteomics-grp")
        grp = os.path.join(self.hiveroot, "quobyte", "proteomics-grp", "fasta")
        os.makedirs(grp)
        os.makedirs(os.path.join(vol, "fasta"))
        for p in (grp, os.path.join(vol, "fasta")):
            open(os.path.join(p, "human.fasta"), "w").close()
        self.set_mount(f"//someone@some-server/proteomics-grp on {vol} (smbfs, nodev)\n")
        j = json.loads(self.run_script(HIVE_PATH, os.path.join(vol, "fasta")).stdout)
        self.assertEqual(j["candidates"], ["/quobyte/proteomics-grp/fasta"])
        self.assertIs(j["verified"], True, j["how"])

    def test_linux_cifs_and_nfs_mounts(self):
        cifs = self._mk("mnt", "prot")
        nfs = self._mk("mnt", "bio")
        self.set_mount(f"//128.120.208.24/proteomics on {cifs} type cifs (rw,relatime)\n"
                       f"10.66.43.50:/flinders/export/group/bioinfo on {nfs} type nfs4 (rw,relatime)\n")
        j = json.loads(self.run_script(HIVE_PATH, os.path.join(cifs, "Data")).stdout)
        self.assertEqual((j["share"], j["rest"]), ("proteomics", "Data"))
        j = json.loads(self.run_script(HIVE_PATH, os.path.join(nfs, "x")).stdout)
        self.assertEqual((j["server"], j["share"]), ("10.66.43.50", "bioinfo"))
        self.assertEqual(j["candidates"], ["/nfs/lssc0/flinders/bioinfo/x", "/quobyte/bioinfo/x"])

    def test_an_ordinary_local_folder_makes_no_ssh_call(self):
        local = self._mk("scripts")
        r = self.run_script(HIVE_PATH, local)
        self.assertEqual(r.returncode, 3, r.stdout + r.stderr)
        self.assertIs(json.loads(r.stdout)["verified"], False)
        self.assertEqual(self.hive_calls(), 0)


class PutGuardTests(Harness):
    def setUp(self):
        super().setUp()
        self.vol = self._mk("Volumes", "proteomics")
        shutil.rmtree(self.vol)
        shutil.copytree(self.shareroot, self.vol)
        self.set_mount(f"//u@128.120.208.24/proteomics on {self.vol} (smbfs, nodev)\n")
        self.run_dir = os.path.join(self.vol, self.service, "P1")

    def test_put_refuses_data_already_on_hive_and_names_the_path(self):
        r = self.run_script(HIVE_EXEC, "--put", self.run_dir, "~/proteomics-pipeline/data/")
        self.assertEqual(r.returncode, 4, r.stdout + r.stderr)
        self.assertIn(f"{FLINDERS}/Data/lab/service/P1", r.stderr)
        self.assertEqual(self.calls("rsync") + self.calls("scp"), [])

    def test_force_uploads_anyway(self):
        r = self.run_script(HIVE_EXEC, "--put", self.run_dir, "~/d/", HIVE_PUT_FORCE="1")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(len(self.calls("rsync")), 1)
        self.assertEqual(self.hive_calls(), 0)

    def test_an_unverified_share_path_still_uploads_but_says_so(self):
        open(os.path.join(self.run_dir, "not_on_hive.raw"), "w").close()
        r = self.run_script(HIVE_EXEC, "--put", self.run_dir, "~/d/")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("no HIVE copy was verified", r.stderr)
        self.assertEqual(len(self.calls("rsync")), 1)

    def test_a_plain_local_path_uploads_without_any_ssh(self):
        """The guard runs before EVERY --put; for an ordinary path it must cost no connection
        (HIVE throttles them)."""
        local = self._mk("scripts")
        r = self.run_script(HIVE_EXEC, "--put", local, "~/proteomics-pipeline/")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.hive_calls(), 0)
        self.assertEqual(self.calls("ssh"), [])
        argv = self.calls("rsync")[0]
        self.assertEqual(argv[-2:], [local, f"tester@{HOST}:~/proteomics-pipeline/"])
        ssh_e = argv[argv.index("-e") + 1]
        self.assertIn("StrictHostKeyChecking=accept-new", ssh_e)
        self.assertIn("ControlMaster=auto", ssh_e)

    def test_a_windows_domain_login_is_refused_before_ssh(self):
        for user in ("AD3+gabrig", "AD3\\gabrig"):
            with self.subTest(user=user):
                r = self.run_script(HIVE_EXEC, "hostname", HIVE_USER=user)
                self.assertEqual(r.returncode, 2)
                self.assertIn("'gabrig'", r.stderr)
        self.assertEqual(self.calls("ssh"), [])

    def test_every_ssh_call_accepts_a_new_host_key(self):
        self.run_script(HIVE_EXEC, "hostname")
        argv = self.calls("ssh")[0]
        self.assertIn("StrictHostKeyChecking=accept-new", argv)
        self.assertEqual(argv[-1], "bash -l -c hostname")


class ScpFallbackTests(Harness):
    SKIP_SHIMS = ("rsync",)          # Git Bash on Windows: no rsync at all

    def setUp(self):
        super().setUp()
        self.src = self._mk("scripts")
        open(os.path.join(self.src, "run_search.py"), "w").close()

    def test_put_a_folder_creates_the_destination_and_lands_inside_it(self):
        r = self.run_script(HIVE_EXEC, "--put", self.src, "~/proteomics-pipeline/")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.calls("ssh")[0][-1], "mkdir -p -- proteomics-pipeline/")
        argv = self.calls("scp")[0]
        self.assertEqual(argv[-2:], [self.src, f"tester@{HOST}:proteomics-pipeline/"])
        for opt in ("IdentitiesOnly=yes", "StrictHostKeyChecking=accept-new", "ControlMaster=auto"):
            self.assertIn(opt, argv)
        self.assertEqual(argv[argv.index("-i") + 1], self.key)
        self.assertIn("-r", argv)

    def test_put_a_folder_to_a_destination_without_a_slash_still_lands_inside_it(self):
        self.run_script(HIVE_EXEC, "--put", self.src, "~/proteomics-pipeline")
        self.assertEqual(self.calls("ssh")[0][-1], "mkdir -p -- proteomics-pipeline")
        self.assertEqual(self.calls("scp")[0][-1], f"tester@{HOST}:proteomics-pipeline/")

    def test_on_windows_there_is_no_controlmaster(self):
        r = self.run_script(HIVE_EXEC, "--put", self.src, "~/x/", FAKE_UNAME="MINGW64_NT-10.0-22631")
        self.assertEqual(r.returncode, 0, r.stderr)
        argv = self.calls("scp")[0]
        self.assertNotIn("ControlMaster=auto", argv)
        self.assertIn("StrictHostKeyChecking=accept-new", argv)

    def test_a_trailing_slash_source_is_refused_not_reinterpreted(self):
        """rsync: dir/ = the CONTENTS. scp copies the folder either way."""
        r = self.run_script(HIVE_EXEC, "--put", self.src + "/", "~/x/")
        self.assertEqual(r.returncode, 2)
        self.assertIn("CONTENTS", r.stderr)
        self.assertEqual(self.calls("scp"), [])

    def test_put_a_file_to_a_new_name_needs_no_mkdir(self):
        f = os.path.join(self.src, "run_search.py")
        self.run_script(HIVE_EXEC, "--put", f, "~/proteomics-pipeline/job.sh")
        self.assertEqual(self.calls("ssh"), [])
        self.assertEqual(self.calls("scp")[0][-1], f"tester@{HOST}:proteomics-pipeline/job.sh")

    def test_get_into_a_folder(self):
        dst = os.path.join(self.d, "session", "output") + "/"
        r = self.run_script(HIVE_EXEC, "--get", "~/proteomics-pipeline/out", dst)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(os.path.isdir(dst))
        self.assertEqual(self.calls("scp")[0][-2:], [f"tester@{HOST}:proteomics-pipeline/out", dst])

    def test_get_one_file_to_a_new_name_leaves_a_file(self):
        """rsync names a single file <dest>; scp is pointed at a fresh <dest>/ and the lone
        file is unwrapped, so the result is the same."""
        dst = os.path.join(self.d, "report.parquet")
        r = self.run_script(HIVE_EXEC, "--get", "~/out/report.parquet", dst,
                            FAKE_SCP_CREATE="report.parquet")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(os.path.isfile(dst))


class CheckAccessTests(Harness):
    def check(self, user="tester", **env):
        r = self.run_script(CHECK_ACCESS, user, self.key, **env)
        self.assertEqual(r.returncode, 0, r.stderr)
        return json.loads(r.stdout)

    def test_the_json_parses_and_has_the_new_fields(self):
        j = self.check(FAKE_KNOWN_FP=PUBLISHED_ED25519)
        self.assertNotIn("proteomics_grp_access", j)            # renamed: it is a LOCAL check
        for k in ("local_proteomics_grp_access", "core_member", "hive_ssh_error",
                  "hive_host_key_known", "hive_host_key_fingerprints", "local_python3"):
            self.assertIn(k, j)
        self.assertEqual(j["hive_ssh"], "ok")
        self.assertIsNone(j["hive_ssh_error"])
        self.assertIs(j["core_member"], True)
        self.assertIs(j["hive_host_key_known"], True)
        self.assertEqual(self.calls("ssh-keyscan"), [])         # known: no extra connection
        self.assertEqual(j["recommended_mode"], "hive_remote")

    def test_group_membership_needs_a_listing_not_ls_d(self):
        """/quobyte is 777 and proteomics-grp is 2770: `ls -d` succeeds for every account."""
        j = self.check(FAKE_SSH="old_ls_d")
        self.assertIs(j["hive_ssh_has_proteomics_grp"], False)
        self.assertEqual(j["core_member"], j["local_proteomics_grp_access"])
        probe = self.calls("ssh")[0][-1]
        self.assertIn("ls /quobyte/proteomics-grp", probe)
        self.assertNotIn("ls -d", probe)

    def test_an_unknown_host_key_is_scanned_and_accepted_when_it_matches(self):
        j = self.check(FAKE_SCAN_FP=PUBLISHED_ED25519)
        self.assertIs(j["hive_host_key_known"], False)
        self.assertEqual(j["hive_host_key_fingerprints"], [f"ED25519 {PUBLISHED_ED25519}"])
        self.assertIs(j["hive_host_key_matches_published"], True)
        self.assertIn("StrictHostKeyChecking=accept-new", self.calls("ssh")[0])

    def test_an_unknown_key_that_does_not_match_is_not_auto_trusted(self):
        j = self.check(FAKE_SCAN_FP="SHA256:somethingElse", FAKE_SSH="hostkey")
        self.assertIs(j["hive_host_key_matches_published"], False)
        self.assertIn("StrictHostKeyChecking=yes", self.calls("ssh")[0])
        self.assertEqual(j["hive_ssh_error"]["kind"], "host_key")
        self.assertIn("published fingerprint", j["hive_ssh_error"]["detail"])

    def test_failures_are_classified(self):
        for mode, kind, text in (("hostkey", "host_key", "Host key verification failed"),
                                 ("denied", "permission_denied", "Permission denied (publickey)"),
                                 ("throttled", "timeout", "kex_exchange_identification"),
                                 # unclassified: ssh's LAST line, never the pre-auth banner
                                 ("other", "other", "Too many authentication failures")):
            with self.subTest(mode=mode):
                j = self.check(FAKE_SSH=mode, FAKE_KNOWN_FP=PUBLISHED_ED25519)
                self.assertEqual(j["hive_ssh"], "failed")
                self.assertEqual(j["hive_ssh_error"]["kind"], kind)
                self.assertIn(text, j["hive_ssh_error"]["detail"])      # not the banner line

    def test_a_windows_domain_login_is_flagged_and_not_tried(self):
        j = self.check(user="AD3+gabrig")
        self.assertIn("'gabrig'", j["hive_user_warning"])
        self.assertEqual(j["hive_ssh"], "failed")
        self.assertEqual(self.calls("ssh"), [])

    def test_the_windows_store_python_stub_is_reported_and_not_run(self):
        apps = self._mk("WindowsApps")
        marker = os.path.join(self.d, "stub_ran")
        _write_exe(os.path.join(apps, "python3"), f'touch "{marker}"\n')
        j = self.check(PATH=f"{apps}:{self.bin}", FAKE_KNOWN_FP=PUBLISHED_ED25519)
        self.assertIs(j["local_python3"]["usable"], False)
        self.assertIn("WindowsApps", j["local_python3"]["path"])
        self.assertFalse(os.path.exists(marker), "the Store stub was executed")

    def test_a_real_python3_is_usable_and_none_is_reported_as_none(self):
        os.symlink(sys.executable, os.path.join(self.bin, "python3"))
        self.assertIs(self.check(FAKE_KNOWN_FP=PUBLISHED_ED25519)["local_python3"]["usable"], True)
        os.remove(os.path.join(self.bin, "python3"))
        j = self.check(FAKE_KNOWN_FP=PUBLISHED_ED25519)
        self.assertIsNone(j["local_python3"]["path"])
        self.assertIs(j["local_python3"]["usable"], False)


class ScriptModeTests(unittest.TestCase):
    def test_new_scripts_are_executable(self):
        """CI requires mode 100755 on every script with a shebang."""
        self.assertTrue(os.access(HIVE_PATH, os.X_OK))


if __name__ == "__main__":
    unittest.main(verbosity=2)
