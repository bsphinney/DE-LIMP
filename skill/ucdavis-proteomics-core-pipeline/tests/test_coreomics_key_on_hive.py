#!/usr/bin/env python3
"""
The CoreOmics key on HIVE (coreomics_key_to_hive.sh), for staff with no usable Python.

Michelle (msalemi, 2026-10-08, Windows + Git Bash, skill 2.11.2) saved her CoreOmics key with the
skill's own line, on her laptop. The CoreOmics steps (check / identify / fetch / bioshare /
email-draft) were documented as LOCAL, but her laptop's python3 is the Microsoft Store stub, so
nothing could run them there -- and HIVE "deliberately has none". Gabriela reported the same on
2026-10-02. Now the steps run on HIVE, from ~/.coreomics_token, and the agent puts the key there
once:

  * the laptop file is ssh's STANDARD INPUT -- the key is never in an argv (ssh's, bash's, the
    remote command's), never printed, and never copied into a temporary file on the laptop;
  * a user without a laptop key pastes it at a hidden prompt (read -rs) in their own Git Bash,
    and it goes the same way;
  * on HIVE (home folders drwxrwsr-x) it lands in a NEW mode-600 file renamed over the name, so a
    link planted there is replaced, not written through;
  * check_access.sh reports the key on HIVE (true / false / mode-wrong) and where the steps run.

No test reaches HIVE: `ssh` is a recorder on PATH that runs the remote command here with HOME set
to a folder standing in for the HIVE home, and the real hive_exec.sh drives it. Every key here is
FAKE.
"""
import json
import os
import shutil
import stat
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
TO_HIVE = os.path.join(SCRIPTS, "coreomics_key_to_hive.sh")
HIVE_EXEC = os.path.join(SCRIPTS, "hive_exec.sh")

KEY = "c0reFAKE" + "9f8e7d6c5b4a39281706f5e4d3c2b1a0"       # never printed, never in an argv
PASTED = "pastedFAKE" + "0a1b2c3d4e5f6a7b8c9d0e1f2a3b"

# ssh stand-in: one call's argv, one arg per line, then the remote command run HERE against the
# stand-in HIVE home (hive_exec.sh sends `bash -l -c <quoted>`; -l dropped so no login profile
# of the test machine runs). Its standard input is what ssh would forward.
FAKE_SSH = r'''#!/bin/bash
for a in "$@"; do printf '%s\n' "$a"; done >> "$LOGS/ssh.log"; echo "--end--" >> "$LOGS/ssh.log"
[ "${FAKE_SSH_DOWN:-0}" = 1 ] && { echo "ssh: connect to host hive.hpc.ucdavis.edu port 22: Operation timed out" >&2; exit 255; }
cmd="${@: -1}"; cmd="${cmd/#bash -l -c /bash -c }"
HOME="$FAKE_HIVE_HOME" exec bash -c "$cmd"
'''
FAKE_UNAME = '#!/bin/bash\necho "${FAKE_UNAME:-Darwin}"\n'
# Git Bash's cygpath -u for the profile folder: C:\Users\x -> the stand-in laptop profile
FAKE_CYGPATH = '#!/bin/bash\necho "$FAKE_PROFILE"\n'


def _exe(path, body):
    with open(path, "w") as fh:
        fh.write(body)
    os.chmod(path, 0o755)


class Harness(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = os.path.realpath(self._tmp.name)
        mk = lambda *p: (os.makedirs(os.path.join(self.d, *p), exist_ok=True), os.path.join(self.d, *p))[1]
        self.bin, self.logs = mk("bin"), mk("logs")
        self.laptop, self.hive, self.tmpdir = mk("laptop"), mk("hive_home"), mk("laptop_tmp")
        _exe(os.path.join(self.bin, "ssh"), FAKE_SSH)
        _exe(os.path.join(self.bin, "uname"), FAKE_UNAME)
        self.sshkey = os.path.join(self.d, "id_test")
        open(self.sshkey, "w").close()
        self.hive_key = os.path.join(self.hive, ".coreomics_token")

    def tearDown(self):
        self._tmp.cleanup()

    def env(self, **extra):
        e = {"PATH": self.bin + os.pathsep + os.environ.get("PATH", "/usr/bin:/bin"),
             "HOME": self.laptop, "LOGS": self.logs, "FAKE_HIVE_HOME": self.hive,
             "HIVE_USER": "tester", "HIVE_KEY": self.sshkey, "HIVE_SSH_MUX": "0",
             "HIVE_ENV_FILE": os.path.join(self.d, "no.env"), "TMPDIR": self.tmpdir}
        e.update(extra)
        return e

    def run_script(self, *args, stdin=None, **env):
        return subprocess.run(["bash", TO_HIVE, *args], input=stdin, capture_output=True, text=True,
                              env=self.env(**env), cwd=self.d, timeout=60)

    def laptop_key(self, text=KEY + "\n", where=None):
        path = os.path.join(where or self.laptop, ".coreomics_token")
        with open(path, "w") as fh:
            fh.write(text)
        os.chmod(path, 0o600)
        return path

    def ssh_argvs(self):
        p = os.path.join(self.logs, "ssh.log")
        if not os.path.exists(p):
            return []
        out, cur = [], []
        with open(p) as fh:
            for line in fh.read().split("\n"):
                if line == "--end--":
                    out.append(cur)
                    cur = []
                elif line or cur:
                    cur.append(line)
        return out

    def assert_never_exposed(self, r, *secrets):
        for s in secrets:
            for part in (s, s[:12], s[-12:]):
                self.assertNotIn(part, r.stdout + r.stderr, "printed")
                for argv in self.ssh_argvs():
                    self.assertFalse(any(part in a for a in argv), f"in ssh argv: {argv}")
        self.assertEqual(os.listdir(self.tmpdir), [], "a temporary file on the laptop")

    def hive_mode(self):
        return stat.S_IMODE(os.stat(self.hive_key).st_mode)

    def read(self, path):
        with open(path) as fh:
            return fh.read()


class StdinTransferTests(Harness):
    def test_the_laptop_key_goes_over_stdin_to_a_mode_600_file(self):
        """Michelle's case: the key saved with the skill's own line, no Python."""
        src = self.laptop_key()
        r = self.run_script()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        j = json.loads(r.stdout)
        self.assertEqual((j["status"], j["coreomics_key_on_hive"], j["hive_mode"]), ("ok", True, "600"))
        self.assertEqual(self.read(self.hive_key), KEY + "\n")
        self.assertEqual(self.hive_mode(), 0o600)
        (argv,) = self.ssh_argvs()                         # ONE connection
        self.assertEqual(argv[-1][:11], "bash -l -c ")      # hive_exec.sh's own call
        self.assert_never_exposed(r, KEY)
        # the laptop copy is kept by default (save_transcript.py redacts what it can read here)
        self.assertEqual((j["laptop_copy"], j["sent_from"]), ("kept", src))
        self.assertTrue(os.path.isfile(src))
        self.assertIn("core_submission.py check --json", j["next"])
        self.assertEqual(sorted(os.listdir(self.hive)), [".coreomics_token"], "no stray file on HIVE")

    def test_remove_local_deletes_the_laptop_copy_only_after_hive_has_it(self):
        src = self.laptop_key()
        r = self.run_script("--remove-local", FAKE_SSH_DOWN="1")
        self.assertEqual(r.returncode, 3, r.stdout + r.stderr)
        self.assertTrue(os.path.isfile(src), "HIVE unreachable: the only copy stays")
        self.assertIn("not saved on HIVE", json.loads(r.stdout)["say"])
        r = self.run_script("--remove-local")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(json.loads(r.stdout)["laptop_copy"], "removed")
        self.assertFalse(os.path.exists(src))
        self.assertEqual(self.read(self.hive_key), KEY + "\n")

    def test_git_bash_reads_the_profile_folder_first(self):
        """On Windows the skill's line saves to $USERPROFILE (Python's ~), which Git Bash's $HOME
        need not be (a domain account)."""
        profile = os.path.join(self.d, "Users", "msalemi")
        os.makedirs(profile)
        self.laptop_key(where=profile)
        self.laptop_key(text="an-older-key-in-git-bash-home\n")
        _exe(os.path.join(self.bin, "cygpath"), FAKE_CYGPATH)
        r = self.run_script(FAKE_UNAME="MINGW64_NT-10.0-22631", USERPROFILE="C:\\Users\\msalemi",
                            FAKE_PROFILE=profile)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(json.loads(r.stdout)["sent_from"], os.path.join(profile, ".coreomics_token"))
        self.assertEqual(self.read(self.hive_key), KEY + "\n")
        self.assert_never_exposed(r, KEY)

    def test_a_readable_key_already_on_hive_is_replaced_by_a_private_file(self):
        with open(self.hive_key, "w") as fh:
            fh.write("old\n")
        os.chmod(self.hive_key, 0o644)
        self.laptop_key()
        r = self.run_script()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.hive_mode(), 0o600)
        self.assertEqual(self.read(self.hive_key), KEY + "\n")

    def test_a_link_planted_at_the_name_is_replaced_not_written_through(self):
        """HIVE home folders are group-writable: a group member can leave a link at
        ~/.coreomics_token to a file of theirs. The key must not go into that file."""
        theirs = os.path.join(self.d, "someone_elses_file")
        with open(theirs, "w") as fh:
            fh.write("")
        os.chmod(theirs, 0o666)
        os.symlink(theirs, self.hive_key)
        self.laptop_key()
        r = self.run_script()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertFalse(os.path.islink(self.hive_key))
        self.assertEqual(self.read(theirs), "", "the key went through the link")
        self.assertEqual(self.read(self.hive_key), KEY + "\n")

    def test_no_laptop_key_sends_nothing_and_names_the_paste_route(self):
        r = self.run_script()
        self.assertEqual(r.returncode, 2, r.stderr)
        j = json.loads(r.stdout)
        self.assertEqual(j["status"], "no_key_here")
        self.assertIn(os.path.join(self.laptop, ".coreomics_token"), j["say"])
        self.assertIn("coreomics_key_to_hive.sh' --paste", j["fix"])
        self.assertIn("own Git Bash window", j["fix"])
        self.assertEqual(self.ssh_argvs(), [], "no connection without a key")

    def test_an_empty_or_huge_file_is_not_sent(self):
        self.laptop_key(text="")
        r = self.run_script()
        self.assertEqual((r.returncode, json.loads(r.stdout)["status"]), (2, "no_key_here"))
        self.laptop_key(text="x" * 5000)
        r = self.run_script()
        self.assertEqual((r.returncode, json.loads(r.stdout)["status"]), (2, "not_a_key"))
        self.assertEqual(self.ssh_argvs(), [])
        self.assertFalse(os.path.exists(self.hive_key))

    def test_the_remote_command_carries_no_secret_and_writes_through_umask_077(self):
        """The command string is the same for every user and every key: it can be logged."""
        r = subprocess.run(["bash", "-c", f'. "{TO_HIVE}"; printf "%s" "$COREOMICS_KEY_STORE"'],
                           capture_output=True, text=True, timeout=30)
        store = r.stdout
        self.assertTrue(store.startswith("umask 077; "))
        self.assertIn("mktemp", store)
        self.assertIn('mv -f -- "$t" "$HOME/.coreomics_token"', store)
        self.assertIn("cat >", store)
        self.assertNotIn("'", r.stdout.split("COREOMICS_KEY_STORE=ok")[1], "probe must splice into '...'")


class PasteTests(Harness):
    def test_a_pasted_key_goes_only_to_hive(self):
        r = self.run_script("--paste", stdin=PASTED + "\n")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        j = json.loads(r.stdout)
        self.assertEqual((j["sent_from"], j["laptop_copy"]), ("the hidden prompt", "none"))
        self.assertEqual(self.read(self.hive_key), PASTED + "\n")
        self.assertEqual(self.hive_mode(), 0o600)
        self.assertEqual(os.listdir(self.laptop), [], "nothing saved on the laptop")
        self.assert_never_exposed(r, PASTED)

    def test_surrounding_blanks_are_dropped(self):
        r = self.run_script("--paste", stdin="  " + PASTED + " \r\n")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.read(self.hive_key), PASTED + "\n")

    def test_nothing_pasted_leaves_the_key_on_hive_alone(self):
        """Ctrl-D, an empty Enter, or the agent running it with no terminal: nothing is sent."""
        with open(self.hive_key, "w") as fh:
            fh.write(KEY + "\n")
        os.chmod(self.hive_key, 0o600)
        for stdin in ("", "\n", "   \n"):
            with self.subTest(stdin=repr(stdin)):
                r = self.run_script("--paste", stdin=stdin)
                self.assertEqual(r.returncode, 2, r.stderr)
                self.assertEqual(json.loads(r.stdout)["status"], "nothing_pasted")
        self.assertEqual(self.ssh_argvs(), [])
        self.assertEqual(self.read(self.hive_key), KEY + "\n")

    def test_remove_local_is_refused_with_paste(self):
        r = self.run_script("--paste", "--remove-local", stdin=PASTED + "\n")
        self.assertEqual(r.returncode, 2)
        self.assertEqual(self.ssh_argvs(), [])


class StatusAndModeTests(Harness):
    """The mode check, from the laptop: ~/.coreomics_token on HIVE must be owner-only."""

    def put(self, mode, text=KEY + "\n"):
        with open(self.hive_key, "w") as fh:
            fh.write(text)
        os.chmod(self.hive_key, mode)

    def status(self):
        r = self.run_script("--status")
        self.assert_never_exposed(r, KEY)
        return r.returncode, json.loads(r.stdout)

    def test_owner_only_is_true(self):
        for mode in (0o600, 0o400):
            with self.subTest(mode=oct(mode)):
                self.put(mode)
                rc, j = self.status()
                self.assertEqual((rc, j["coreomics_key_on_hive"], j["hive_mode"]), (0, True, f"{mode:o}"))
                os.chmod(self.hive_key, 0o600)

    def test_any_group_or_other_bit_is_mode_wrong_with_the_chmod(self):
        for mode in (0o644, 0o640, 0o604, 0o660, 0o666):
            with self.subTest(mode=oct(mode)):
                self.put(mode)
                rc, j = self.status()
                self.assertEqual((rc, j["coreomics_key_on_hive"], j["hive_mode"]), (2, "mode-wrong", f"{mode:o}"))
                self.assertIn("other HIVE accounts can read it", j["say"])
                self.assertIn("hive_exec.sh 'chmod 600 ~/.coreomics_token'", j["fix"])

    def test_absent_empty_and_not_a_file_are_false(self):
        rc, j = self.status()
        self.assertEqual((rc, j["coreomics_key_on_hive"], j["status"]), (2, False, "absent"))
        self.put(0o600, text="")
        self.assertIs(self.status()[1]["coreomics_key_on_hive"], False)
        os.remove(self.hive_key)
        os.makedirs(self.hive_key)
        self.assertIs(self.status()[1]["coreomics_key_on_hive"], False)

    def test_unreachable_is_exit_3_not_absent(self):
        r = self.run_script("--status", FAKE_SSH_DOWN="1")
        self.assertEqual(r.returncode, 3)
        self.assertIsNone(json.loads(r.stdout)["coreomics_key_on_hive"])


class RouteTests(unittest.TestCase):
    """Where the CoreOmics steps run (coreomics_route, which check_access.sh reports)."""

    def route(self, mode, py, here, hive):
        r = subprocess.run(["bash", "-c", f'. "{TO_HIVE}"; coreomics_route "$@"', "_",
                            mode, py, here, hive], capture_output=True, text=True, timeout=30)
        self.assertEqual(r.returncode, 0, r.stderr)
        return r.stdout.strip()

    def test_the_table(self):
        cases = [
            # mode          python  key here  key on HIVE   -> runs on
            ("hive_remote", "false", "true", "false", "hive"),        # Michelle, before setup
            ("hive_remote", "false", "true", "true", "hive"),         # ... and after
            ("hive_remote", "false", "false", "null", "hive"),
            ("hive_remote", "false", "false", "mode-wrong", "hive"),
            ("hive_remote", "true", "true", "true", "this_computer"),  # a Mac with Python: unchanged
            ("hive_remote", "true", "true", "false", "this_computer"),
            ("hive_remote", "true", "false", "true", "hive"),          # key moved to HIVE only
            ("hive_remote", "true", "false", "mode-wrong", "hive"),
            ("hive_remote", "true", "false", "false", "this_computer"),
            ("local", "false", "false", "null", "this_computer"),      # no HIVE: never HIVE
            ("local", "true", "true", "null", "this_computer"),
            ("hive_local", "true", "true", "true", "this_computer"),   # Claude Code on HIVE itself
        ]
        for mode, py, here, hive, want in cases:
            with self.subTest(mode=mode, py=py, here=here, hive=hive):
                self.assertEqual(self.route(mode, py, here, hive), want)

    def test_sourcing_runs_nothing(self):
        r = subprocess.run(["bash", "-c", f'. "{TO_HIVE}"; echo sourced'], capture_output=True,
                           text=True, timeout=30, input="")
        self.assertEqual((r.returncode, r.stdout, r.stderr), (0, "sourced\n", ""))

    def test_the_probe_answers_without_reading_the_key(self):
        with tempfile.TemporaryDirectory() as home:
            def probe():
                r = subprocess.run(["bash", "-c", f'. "{TO_HIVE}"; bash -c "$COREOMICS_KEY_PROBE"'],
                                   capture_output=True, text=True, timeout=30, env=dict(os.environ, HOME=home))
                self.assertNotIn(KEY, r.stdout + r.stderr)
                return r.stdout.strip()
            self.assertEqual(probe(), "COREOMICS_KEY=absent")
            f = os.path.join(home, ".coreomics_token")
            with open(f, "w") as fh:
                fh.write(KEY)
            os.chmod(f, 0o640)
            self.assertEqual(probe(), "COREOMICS_KEY=mode:640")
            os.chmod(f, 0o600)
            self.assertEqual(probe(), "COREOMICS_KEY=mode:600")


class ScriptModeTests(unittest.TestCase):
    def test_it_is_executable_and_parses(self):
        """CI requires 100755 on every script with a shebang; SKILL.md runs it with bash."""
        self.assertTrue(os.access(TO_HIVE, os.X_OK))
        self.assertEqual(subprocess.run(["bash", "-n", TO_HIVE]).returncode, 0)
        r = subprocess.run(["bash", TO_HIVE, "--help"], capture_output=True, text=True, timeout=30)
        self.assertEqual(r.returncode, 0)
        self.assertIn("--paste", r.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
