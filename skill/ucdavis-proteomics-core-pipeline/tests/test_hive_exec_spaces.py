#!/usr/bin/env python3
"""
hive_exec.sh --get/--put with a remote path that has a space (2.14.1 I16/I33): rsync and
scp re-split the path on the remote side, so "UC Davis Lab/logs" became a stray file "UC".
Such a path is streamed as tar over ssh. HIVE is a directory tree: the fake ssh runs the
remote command locally. Plain paths must keep their rsync command.
"""
import os
import shutil
import stat
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
HIVE_EXEC = os.path.join(os.path.dirname(HERE), "scripts", "hive_exec.sh")
TOOLS = ("tar", "mkdir", "cp", "rm", "mktemp", "dirname", "basename", "cat", "ls",
         "tr", "uname", "mv", "rmdir", "sed", "grep", "head", "awk", "sort", "wc", "chmod")
SSH = '#!/bin/bash\nfor last; do :; done\nexec bash -c "$last"\n'
RSYNC = '#!/bin/bash\nfor a in "$@"; do printf "%s\\n" "$a"; done >> "$LOGS/rsync.log"\necho --end-- >> "$LOGS/rsync.log"\n'


class Spaces(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = os.path.realpath(self._t.name)
        self.bin = self._mk("bin")
        self.logs = self._mk("logs")
        self.remote = self._mk("q", "off_campus", "UC Davis Lab", "logs", "conversation")
        self.local = self._mk("local")
        for t in TOOLS:
            real = shutil.which(t)
            if not real:
                self.skipTest(f"{t} not installed")
            os.symlink(real, os.path.join(self.bin, t))
        os.symlink("/bin/bash", os.path.join(self.bin, "bash"))   # macOS bash 3.2
        self._exe("ssh", SSH)
        self.key = os.path.join(self.d, "k")
        open(self.key, "w").close()

    def tearDown(self):
        self._t.cleanup()

    def _mk(self, *p):
        q = os.path.join(self.d, *p)
        os.makedirs(q, exist_ok=True)
        return q

    def _exe(self, name, body):
        p = os.path.join(self.bin, name)
        with open(p, "w") as fh:
            fh.write(body)
        os.chmod(p, os.stat(p).st_mode | stat.S_IXUSR)

    def run_exec(self, *args, cwd=None):
        env = {"PATH": self.bin, "HOME": self._mk("home"), "LOGS": self.logs, "HIVE_USER": "tester",
               "HIVE_KEY": self.key, "HIVE_ENV_FILE": os.path.join(self.d, "no.env"),
               "TMPDIR": self._mk("tmp"), "HIVE_PUT_FORCE": "1", "HIVE_SSH_MUX": "0"}
        return subprocess.run(["/bin/bash", HIVE_EXEC, *args], capture_output=True, text=True,
                              env=env, cwd=cwd or self.d, timeout=60)

    def write(self, path, text):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)

    def test_get_from_a_folder_with_a_space(self):
        self.write(os.path.join(self.remote, "a.jsonl"), "A")
        self.write(os.path.join(self.remote, "sub", "b.jsonl"), "B")
        dest = self._mk("tmpdest")
        r = self.run_exec("--get", self.remote, dest)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(open(os.path.join(dest, "conversation", "a.jsonl")).read(), "A")
        self.assertEqual(open(os.path.join(dest, "conversation", "sub", "b.jsonl")).read(), "B")

    def test_get_trailing_slash_copies_the_contents(self):
        self.write(os.path.join(self.remote, "a.jsonl"), "A")
        dest = os.path.join(self.d, "newdest")
        r = self.run_exec("--get", self.remote + "/", dest)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(os.path.exists(os.path.join(dest, "a.jsonl")))

    def test_put_into_a_folder_with_a_space(self):
        src = os.path.join(self.local, "x.jsonl")
        self.write(src, "X")
        dest = os.path.join(self.d, "q", "off_campus", "UC Davis Lab", "logs", "conversation") + "/"
        r = self.run_exec("--put", src, dest)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(open(os.path.join(self.remote, "x.jsonl")).read(), "X")
        self.assertFalse(os.path.exists(os.path.join(self.d, "q", "off_campus", "UC")))

    def test_put_a_folder_lands_inside_a_new_dest(self):
        self.write(os.path.join(self.local, "pkg", "s.sh"), "S")
        dest = os.path.join(self.d, "q", "new dir")
        r = self.run_exec("--put", os.path.join(self.local, "pkg"), dest + "/")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(open(os.path.join(dest, "pkg", "s.sh")).read(), "S")

    def test_path_with_a_quote_stays_one_word(self):
        odd = self._mk("q", "it's here")
        self.write(os.path.join(odd, "f.txt"), "F")
        dest = self._mk("d2")
        r = self.run_exec("--get", odd, dest)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(os.path.exists(os.path.join(dest, "it's here", "f.txt")))

    def test_single_file_to_a_missing_dest_is_a_file(self):
        f = os.path.join(self.remote, "one.txt")
        self.write(f, "1")
        dest = os.path.join(self.d, "got.txt")
        r = self.run_exec("--get", f, dest)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(os.path.isfile(dest))
        # put: a file to a missing remote name becomes that file
        rdest = os.path.join(self.remote, "renamed.txt")
        r = self.run_exec("--put", dest, rdest)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(open(rdest).read(), "1")

    def test_get_of_a_missing_remote_fails(self):
        r = self.run_exec("--get", os.path.join(self.remote, "nope"), self._mk("d3"))
        self.assertNotEqual(r.returncode, 0)

    def test_a_login_banner_on_stdout_does_not_corrupt_the_stream(self):
        self._exe("ssh", SSH.replace("exec bash", 'case "$last" in "bash -l"*) echo BANNER ;; esac\nexec bash'))
        self.write(os.path.join(self.remote, "a.jsonl"), "A")
        dest = self._mk("d4")
        r = self.run_exec("--get", self.remote, dest)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(open(os.path.join(dest, "conversation", "a.jsonl")).read(), "A")

    def test_names_starting_with_a_dash(self):
        dashdir = self._mk("q", "-dash dir")
        self.write(os.path.join(dashdir, "f.txt"), "F")
        dest = self._mk("d5")
        r = self.run_exec("--get", dashdir, dest)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(os.path.exists(os.path.join(dest, "-dash dir", "f.txt")))
        self.write(os.path.join(self.local, "-f x"), "X")
        rdest = os.path.join(self.remote, "dd") + "/"
        r = self.run_exec("--put", "-f x", rdest, cwd=self.local)
        self.assertEqual(r.returncode, 0, r.stderr)
        r = self.run_exec("--put", "./-f x", rdest, cwd=self.local)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(open(os.path.join(self.remote, "dd", "-f x")).read(), "X")

    def test_put_of_an_executable_file_to_a_new_name_stays_executable(self):
        src = os.path.join(self.local, "run.sh")
        self.write(src, "#!/bin/sh\n")
        os.chmod(src, 0o755)
        rdest = os.path.join(self.remote, "run copy.sh")
        r = self.run_exec("--put", src, rdest)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(os.access(rdest, os.X_OK))

    def test_get_into_an_existing_dest_keeps_mtimes_and_leaves_no_stage(self):
        f = os.path.join(self.remote, "a.jsonl")
        self.write(f, "A")
        os.utime(f, (1000000000, 1000000000))
        dest = self._mk("d6")
        r = self.run_exec("--get", self.remote, dest)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(int(os.stat(os.path.join(dest, "conversation", "a.jsonl")).st_mtime), 1000000000)
        self.assertEqual(os.listdir(dest), ["conversation"])
        self.assertEqual([x for x in os.listdir(self.d) if x.startswith(".hive_get")], [])

    def test_a_failed_get_into_an_existing_dest_changes_nothing(self):
        dest = self._mk("d7")
        r = self.run_exec("--get", os.path.join(self.remote, "nope"), dest + "/")
        self.assertNotEqual(r.returncode, 0)
        self.assertEqual(os.listdir(dest), [])

    def test_plain_paths_still_use_rsync(self):
        self._exe("rsync", RSYNC)
        plain = os.path.join(self.d, "plainremote")
        r = self.run_exec("--get", plain, self.local)
        self.assertEqual(r.returncode, 0, r.stderr)
        log = open(os.path.join(self.logs, "rsync.log")).read().split("\n")
        self.assertIn("-rlt", log)
        self.assertTrue(any(a.endswith("@hive.hpc.ucdavis.edu:" + plain) for a in log), log)


if __name__ == "__main__":
    unittest.main()
