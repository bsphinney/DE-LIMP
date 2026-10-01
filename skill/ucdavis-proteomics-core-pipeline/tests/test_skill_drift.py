#!/usr/bin/env python3
"""
Skill drift: HIVE must run the skill this computer runs, and staff must hear when a newer
release is out (skill_version.sh --check-hive, SKILL.md step 0).

On 2026-09-29 two Core staff were still on 2.6.0 -- on their laptops AND in ~/proteomics-pipeline
on HIVE -- while 2.8.0 was out, so none of the 2.7/2.8 fixes had ever run for them. One HIVE copy
had no .claude-plugin/ at all (its sdrf.tsv said "v0.0.0"). Nothing compared the copies and
nothing announced the release.

No test reaches HIVE: HIVE_EXEC is a stand-in that runs the command here with HOME set to a
folder standing in for the HIVE home (a symlink, so its `pwd -P` form differs), `squeue` is a
stand-in (SKILL_SQUEUE) that prints the batch-script paths listed in $FAKE_SQUEUE_LIST -- or
fails like a timed-out slurmctld with FAKE_SQUEUE_FAIL=1 -- the release folder is a temporary
one, and a separate case drives the real hive_exec.sh through a fake `ssh`.
"""
import json
import os
import pwd
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, HERE)
from job_env import job_env     # noqa: E402

PUT = "bash scripts/hive_exec.sh --put-skill"
UPDATE = "claude plugin update ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core"


def sh_exe(path, body):
    with open(path, "w") as fh:
        fh.write("#!/usr/bin/env bash\n" + body)
    os.chmod(path, 0o755)
    return path


class Harness(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.laptop = os.path.join(self.d, "laptop")       # the installed skill on this computer
        os.makedirs(os.path.join(self.laptop, "scripts"))
        for f in ("skill_version.sh", "hive_exec.sh", "hive_path.sh", "hive_shares.tsv"):
            shutil.copy(os.path.join(SCRIPTS, f), os.path.join(self.laptop, "scripts"))
        # the admins: the shipped list plus whoever runs the suite -- `id -un` on the stand-in
        # HIVE -- so the publish tests can publish (PublishRelease.test_only_an_admin_... not)
        self.me = pwd.getpwuid(os.getuid()).pw_name
        self.admins([self.me])
        # $HOME on "HIVE", reached through a symlink (as /home -> /quobyte/home can be)
        self.hive_home_real = os.path.realpath(os.path.join(self.d, "hive_home_real"))
        os.makedirs(self.hive_home_real)
        self.hive_home = os.path.join(self.d, "hive_home")
        os.symlink(self.hive_home_real, self.hive_home)
        self.hive_root = os.path.join(self.hive_home, "proteomics-pipeline")
        self.grp = os.path.join(self.d, "grp")             # /quobyte/proteomics-grp
        self.release_dir = os.path.join(self.grp, "skill_release")
        os.makedirs(self.grp)
        self.key = os.path.join(self.d, "id_test")
        open(self.key, "w").close()
        self.calls = os.path.join(self.d, "hive_calls.log")
        # hive_exec.sh stand-in: one line per call; FAKE_HIVE_DOWN=1 fails like a dead link, in
        # colour, as some terminals' ssh does
        self.fake = sh_exe(os.path.join(self.d, "fake_hive_exec.sh"),
                           'echo call >> "$FAKE_CALLS"\n'
                           'if [ "${FAKE_HIVE_DOWN:-0}" = 1 ]; then\n'
                           '  printf "\\033[31mssh: connect to host hive.hpc.ucdavis.edu port 22: '
                           'Operation timed out\\033[0m\\a\\n" >&2\n'
                           '  exit 255\nfi\n'
                           'if [ "${FAKE_HIVE_DOWN:-0}" = latin1 ]; then\n'
                           '  printf "ssh: Verbindung zu \\374berpr\\374ft \\377\\376 fehlgeschlagen\\n" >&2\n'
                           '  exit 255\nfi\n'
                           'HOME="$FAKE_HIVE_HOME" exec bash -c "$1"\n')
        # who owns a path on "HIVE": the release file FAKE_FILE_OWNER, anything else (the
        # folder) FAKE_DIR_OWNER -- an admin unless a test says otherwise
        self.owner = sh_exe(os.path.join(self.d, "fake_owner"),
                            'case "$1" in */CURRENT_VERSION) echo "${FAKE_FILE_OWNER-brettsp}" ;;\n'
                            '  *) echo "${FAKE_DIR_OWNER-brettsp}" ;; esac\n')
        # squeue -h -u <user> -o %o: one batch-script path (or command) per job
        self.squeue = sh_exe(os.path.join(self.d, "fake_squeue"),
                             'if [ "${FAKE_SQUEUE_FAIL:-0}" = 1 ]; then\n'
                             '  echo "squeue: error: slurm_load_jobs error: Socket timed out on '
                             'send/recv operation" >&2; exit 1\nfi\n'
                             'case "$*" in *"-o %o"*) ;; *) echo "not asked for %o: $*" >&2; exit 2 ;; esac\n'
                             '[ -n "${FAKE_SQUEUE_LIST:-}" ] && cat "$FAKE_SQUEUE_LIST"\n'
                             'exit 0\n')

    def admins(self, extra=(), text=None):
        """This computer's scripts/core_admins.txt: the shipped one plus `extra` (or `text`)."""
        if text is None:
            with open(os.path.join(SCRIPTS, "core_admins.txt"), encoding="utf-8") as fh:
                text = fh.read() + "".join(f"{u}\n" for u in extra)
        with open(os.path.join(self.laptop, "scripts", "core_admins.txt"), "w") as fh:
            fh.write(text)

    def tearDown(self):
        for root, dirs, _ in os.walk(self.d):
            for x in dirs:
                os.chmod(os.path.join(root, x), 0o755)
        self._tmp.cleanup()

    # ---------------------------------------------------------------- fixtures
    def plugin(self, root, version):
        """<root>/.claude-plugin/plugin.json with this version (None: no .claude-plugin/)."""
        if version is None:
            shutil.rmtree(os.path.join(root, ".claude-plugin"), ignore_errors=True)
            return
        os.makedirs(os.path.join(root, ".claude-plugin"), exist_ok=True)
        with open(os.path.join(root, ".claude-plugin", "plugin.json"), "w") as fh:
            json.dump({"name": "ucdavis-proteomics-core-pipeline", "version": version}, fh)

    def hive_copy(self, version, scripts=True):
        """~/proteomics-pipeline on HIVE as --put-skill leaves it: this computer's scripts/ and a
        plugin.json with `version` (None: no .claude-plugin/, as msalemi's copy had)."""
        if scripts:
            shutil.rmtree(os.path.join(self.hive_root, "scripts"), ignore_errors=True)
            shutil.copytree(os.path.join(self.laptop, "scripts"),
                            os.path.join(self.hive_root, "scripts"))
        self.plugin(self.hive_root, version)

    def release(self, text):
        os.makedirs(self.release_dir, exist_ok=True)
        with open(os.path.join(self.release_dir, "CURRENT_VERSION"), "w") as fh:
            fh.write(text)

    def env(self, login=True, **extra):
        e = job_env(self.d, HOME=os.path.join(self.d, "laptop_home"),
                    SKILL_RELEASE_DIR=self.release_dir, SKILL_SQUEUE=self.squeue,
                    SKILL_OWNER_CMD=self.owner,
                    FAKE_HIVE_HOME=self.hive_home, FAKE_CALLS=self.calls, TMPDIR=self.d,
                    LC_ALL="C.UTF-8", LANG="C.UTF-8")
        if login:
            e.update(HIVE_USER="gabrig", HIVE_KEY=self.key, HIVE_EXEC=self.fake)
        e.update({k: str(v) for k, v in extra.items()})
        return e

    def run_sv(self, *args, login=True, **extra):
        return subprocess.run(["bash", os.path.join(self.laptop, "scripts", "skill_version.sh"),
                               *args], capture_output=True, text=True, timeout=60,
                              env=self.env(login=login, **extra))

    def check(self, want_rc, *args, **extra):
        r = self.run_sv("--check-hive", *args, **extra)
        self.assertEqual(r.returncode, want_rc, r.stdout + r.stderr)
        return json.loads(r.stdout)

    def n_calls(self):
        if not os.path.exists(self.calls):
            return 0
        with open(self.calls) as fh:
            return len(fh.read().split())

    def jobs(self, *scripts):
        """squeue's %o for the user's jobs: a path per job (made as a batch script naming
        `body`), or a bare command. Returns the FAKE_SQUEUE_LIST file."""
        lines = []
        for i, item in enumerate(scripts):
            if isinstance(item, tuple):             # (body, times): one script, several jobs
                body, times = item
            else:
                body, times = item, 1
            if body is None or body.startswith("CMD:"):
                lines += [(body or "CMD:(null)")[4:]] * times
                continue
            path = os.path.join(self.d, "jobs", f"job{i}.sbatch")
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as fh:
                fh.write("#!/bin/bash\n#SBATCH -p high\n" + body + "\n")
            lines += [path] * times
        lst = os.path.join(self.d, "squeue_list.txt")
        with open(lst, "w") as fh:
            fh.write("".join(l + "\n" for l in lines))
        return lst


class CheckHive(Harness):
    def test_in_step_says_nothing_and_costs_one_call(self):
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.9.0")
        self.release("2.9.0\n# comment\n")
        j = self.check(0, "--mode", "hive_remote")
        self.assertEqual((j["status"], j["local"], j["hive"], j["release"]),
                         ("in_step", "2.9.0", "2.9.0", "2.9.0"))
        self.assertEqual((j["release_state"], j["behind_release"], j["jobs"], j["jobs_total"],
                          j["next"], j["say"]), ("published", False, 0, 0, None, None))
        self.assertEqual(j["hive_path"], "~/proteomics-pipeline")
        self.assertEqual(self.n_calls(), 1)

    def test_a_hive_copy_without_claude_plugin_is_put_again(self):
        """msalemi's HIVE copy: scripts/ but no .claude-plugin/ -- every record said 'unknown'."""
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy(None)
        j = self.check(3)
        self.assertEqual((j["status"], j["hive"], j["next"]), ("hive_missing", None, PUT))
        self.assertIn("unknown", j["say"])
        self.assertIn("Putting this computer's 2.9.0 there", j["say"])

    def test_no_copy_on_hive_at_all(self):
        self.plugin(self.laptop, "2.9.0")
        j = self.check(3)
        self.assertEqual((j["status"], j["next"]), ("hive_missing", PUT))

    def test_a_plugin_json_without_scripts_is_not_a_copy(self):
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.9.0", scripts=False)
        self.assertEqual(self.check(3)["status"], "hive_missing")

    def test_an_older_hive_copy_is_put_again(self):
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.6.0")
        j = self.check(3)
        self.assertEqual((j["status"], j["local"], j["hive"], j["next"]),
                         ("hive_behind", "2.9.0", "2.6.0", PUT))
        self.assertIn("2.6.0", j["say"])

    def test_numbers_compare_as_numbers(self):
        self.plugin(self.laptop, "2.10.0")
        self.hive_copy("2.9.0")
        self.assertEqual(self.check(3)["status"], "hive_behind")    # 2.9.0 < 2.10.0
        self.hive_copy("2.10.0")
        self.plugin(self.laptop, "2.9.0")
        self.assertEqual(self.check(4)["status"], "hive_ahead")

    def test_a_pre_release_compares_on_its_release_number(self):
        """Review: laptop 2.9.0-dev vs HIVE 2.10.0 read as 'different, laptop wins' and put the
        older build over the newer one. The release number decides first."""
        cases = [("2.9.0-dev", "2.10.0", 4, "hive_ahead"),     # never downgrade HIVE
                 ("2.10.0-dev", "2.9.0", 3, "hive_behind"),
                 ("2.10.0", "2.9.0+build7", 3, "hive_behind"),
                 ("2.9.0", "2.9.0-dev", 3, "hive_differs"),    # same release, another build:
                 ("2.9.0-dev", "2.9.0", 3, "hive_differs")]    # this computer's wins
        for here, hive, rc, status in cases:
            with self.subTest(here=here, hive=hive):
                self.plugin(self.laptop, here)
                self.hive_copy(hive)
                self.assertEqual(self.check(rc)["status"], status)

    def test_a_different_build_under_the_same_version_is_put_again(self):
        """Same version string, different files (a test build, a hand edit): the digest says so."""
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.9.0")
        with open(os.path.join(self.hive_root, "scripts", "hive_path.sh"), "a") as fh:
            fh.write("# edited on HIVE\n")
        j = self.check(3)
        self.assertEqual((j["status"], j["next"]), ("hive_differs", PUT))
        self.assertIn("files differ", j["say"])

    def test_a_file_name_with_a_space_is_one_file(self):
        """Review: the digest's file list was split on spaces, so "session 2.py" (an iCloud
        duplicate) read as two missing files on both sides, and HIVE lacking it went unseen."""
        self.plugin(self.laptop, "2.9.0")
        with open(os.path.join(self.laptop, "scripts", "session 2.py"), "w") as fh:
            fh.write("# a duplicate\n")
        self.hive_copy("2.9.0")
        self.assertEqual(self.check(0)["status"], "in_step")
        os.remove(os.path.join(self.hive_root, "scripts", "session 2.py"))
        self.assertEqual(self.check(3)["status"], "hive_differs")

    def test_files_an_older_put_left_on_hive_do_not_count(self):
        """--put-skill never deletes: a script this copy no longer has stays on HIVE forever."""
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.9.0")
        open(os.path.join(self.hive_root, "scripts", "retired_script.py"), "w").close()
        os.makedirs(os.path.join(self.hive_root, "scripts", "__pycache__"))
        self.assertEqual(self.check(0)["status"], "in_step")

    def test_running_jobs_that_use_the_copy_stop_an_automatic_put(self):
        """Job-end hooks run ~/proteomics-pipeline/scripts hours later: never swap them."""
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.6.0")
        lst = self.jobs(f"python3 {self.hive_root}/scripts/notify_slack.py job-end",
                        ("diann --f a.d --out x", 3))
        j = self.check(6, FAKE_SQUEUE_LIST=lst)
        self.assertEqual((j["status"], j["jobs"], j["jobs_total"], j["next"]),
                         ("hive_behind", 1, 4, None))
        self.assertIn("1 of your 4 SLURM job(s)", j["say"])
        self.assertIn(PUT, j["say"])
        self.assertNotIn("Putting this computer", j["say"])

    def test_only_jobs_whose_scripts_use_the_copy_count(self):
        """Review: 228 of Brett's jobs, 1 of which used the copy -- counting all of them
        stopped every session with a message false for 227. Every spelling of the copy's
        scripts/ counts; its lookalikes do not; an unreadable script counts (fail closed)."""
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.6.0")
        real = os.path.realpath(self.hive_root)
        self.assertNotEqual(real, self.hive_root)          # the pwd -P form is a different text
        uses = [f"python3 {self.hive_root}/scripts/record_run.py search-done",
                f"python3 {real}/scripts/notify_slack.py",
                "python3 ~/proteomics-pipeline/scripts/fran_deposit.py stage",
                "bash $HOME/proteomics-pipeline/scripts/watch_run.sh",
                'python3 "${HOME}/proteomics-pipeline/scripts/notify_slack.py"',
                ("python3 ~/proteomics-pipeline/scripts/notify_slack.py", 2)]   # an array: 2 jobs
        not_uses = ["source ~/.proteomics-pipeline/activate.sh",
                    "python3 ~/proteomics-pipeline-280-test/scripts/run_search.py",
                    "cd ~/proteomics-pipeline && ls",
                    "CMD:(null)", "CMD:bash"]              # no script path: an interactive srun
        lst = self.jobs(*uses, *not_uses, "UNREADABLE")
        unreadable = os.path.join(self.d, "jobs", f"job{len(uses) + len(not_uses)}.sbatch")
        os.chmod(unreadable, 0)
        j = self.check(6, FAKE_SQUEUE_LIST=lst)
        self.assertEqual((j["jobs"], j["jobs_total"]), (6 + 1 + 1, 6 + 1 + 5 + 1))
        lst = self.jobs(*not_uses)
        j = self.check(3, FAKE_SQUEUE_LIST=lst)
        self.assertEqual((j["jobs"], j["jobs_total"], j["next"]), (0, 5, PUT))

    def test_jobs_do_not_block_when_nothing_would_be_put(self):
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.9.0")
        lst = self.jobs(f"python3 {self.hive_root}/scripts/notify_slack.py")
        self.assertEqual(self.check(0, FAKE_SQUEUE_LIST=lst)["jobs"], 1)
        self.hive_copy("2.10.0")
        self.assertEqual(self.check(4, FAKE_SQUEUE_LIST=lst)["status"], "hive_ahead")

    def test_a_squeue_that_errors_fails_closed(self):
        """Review: squeue present but failing ("Socket timed out") was read as 'no jobs' and the
        copy was put. It could not count them: exit 6."""
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.6.0")
        j = self.check(6, FAKE_SQUEUE_FAIL=1)
        self.assertEqual((j["status"], j["jobs"], j["jobs_total"], j["next"]),
                         ("hive_behind", None, None, None))
        self.assertIn("Socket timed out", j["say"])
        self.assertIn("could not be counted", j["say"])

    def test_no_squeue_at_all_puts_without_counting(self):
        """Only a HIVE with no squeue is 'unknown and allowed'."""
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.6.0")
        j = self.check(3, SKILL_SQUEUE=os.path.join(self.d, "no_such_squeue"))
        self.assertEqual((j["jobs"], j["jobs_total"], j["next"]), (None, None, PUT))

    def test_a_newer_hive_copy_stops_and_is_never_overwritten_unasked(self):
        self.plugin(self.laptop, "2.6.0")
        self.hive_copy("2.9.0")
        self.release("2.9.0\n")
        j = self.check(4)
        self.assertEqual((j["status"], j["next"], j["behind_release"]), ("hive_ahead", None, True))
        self.assertIn(UPDATE, j["say"])
        self.assertIn("2.9.0", j["say"])

    def test_a_newer_unreleased_hive_build_has_nothing_to_update_to(self):
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.10.0")
        self.release("2.9.0\n")
        j = self.check(4)
        self.assertIs(j["behind_release"], False)
        self.assertIn("nothing to update to", j["say"])
        self.assertIn(PUT, j["say"])
        self.assertNotIn(UPDATE, j["say"])
        j = self.check(4, SKILL_RELEASE_DIR=os.path.join(self.d, "nowhere", "skill_release"))
        self.assertIsNone(j["behind_release"])
        self.assertIn("test build", j["say"])

    def test_a_laptop_behind_the_release_is_warned_with_the_update_that_works(self):
        """gabrig: laptop 2.6.0, HIVE 2.6.0 (in step with each other), release 2.9.0. The update
        is `claude plugin update` -- `/plugin marketplace update` only refreshes the catalogue
        (checked on Claude Code 2.1.285)."""
        self.plugin(self.laptop, "2.6.0")
        self.hive_copy("2.6.0")
        self.release("2.9.0\n")
        j = self.check(0)
        self.assertEqual((j["status"], j["release"], j["behind_release"]),
                         ("in_step", "2.9.0", True))
        self.assertIn("current release is 2.9.0", j["say"])
        self.assertIn(UPDATE, j["say"])
        self.assertIn("Enable auto-update", j["say"])
        self.assertNotIn("marketplace update", j["say"])

    def test_both_at_once_are_both_said(self):
        self.plugin(self.laptop, "2.6.0")
        self.hive_copy(None)
        self.release("2.9.0\n")
        j = self.check(3)
        self.assertEqual((j["status"], j["behind_release"]), ("hive_missing", True))
        self.assertIn("Putting this computer's 2.6.0", j["say"])
        self.assertIn("current release is 2.9.0", j["say"])

    def test_numbers_compare_as_numbers_against_the_release(self):
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.9.0")
        self.release("2.10.0\n")
        self.assertIs(self.check(0)["behind_release"], True)
        self.release("2.8.1\n")
        j = self.check(0)
        self.assertIs(j["behind_release"], False, "ahead of the release (a test build) is fine")
        self.assertIsNone(j["say"])
        self.plugin(self.laptop, "2.9.0-dev")
        self.hive_copy("2.9.0-dev")
        self.release("2.9.0\n")
        self.assertIs(self.check(0)["behind_release"], False, "a build of the release itself")

    def test_the_release_file_is_read_sensibly(self):
        self.plugin(self.laptop, "2.6.0")
        self.hive_copy("2.6.0")
        for text in ("\n\n2.9.0\n", "# written by x\nv2.9.0\n", "  2.9.0  \n# c\n"):
            with self.subTest(text=text):
                self.release(text)
                j = self.check(0)
                self.assertEqual((j["release"], j["release_state"]), ("2.9.0", "published"))

    def test_no_release_published_yet(self):
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.9.0")
        j = self.check(0)
        self.assertEqual((j["release"], j["release_state"], j["behind_release"], j["say"]),
                         (None, "not_published", None, None))

    def test_no_access_to_the_group_folder_is_quiet(self):
        """A HIVE account outside proteomics-grp cannot read the release: nothing is said."""
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.9.0")
        self.release("3.0.0\n")
        os.chmod(self.grp, 0)
        j = self.check(0)
        self.assertEqual((j["release"], j["release_state"], j["say"]), (None, "no_access", None))

    def test_a_garbled_release_file_is_not_a_version(self):
        self.plugin(self.laptop, "2.6.0")
        self.hive_copy("2.6.0")
        self.release("next week\n")
        j = self.check(0)
        self.assertEqual((j["release"], j["release_state"], j["say"]), (None, "unreadable", None))

    def test_hive_unreachable_is_not_fatal_and_its_colours_do_not_break_the_json(self):
        self.plugin(self.laptop, "2.9.0")
        r = self.run_sv("--check-hive", FAKE_HIVE_DOWN=1)
        self.assertEqual(r.returncode, 5, r.stderr)
        j = json.loads(r.stdout)                         # parses: no raw escape characters
        self.assertEqual(j["status"], "unreachable")
        self.assertIn("Operation timed out", j["say"])
        self.assertIn("carrying on", j["say"])
        self.assertFalse(any(ord(c) < 32 for c in r.stdout.rstrip("\n")), repr(r.stdout))
        self.assertNotIn("[31m", j["say"])

    def test_a_non_utf8_ssh_error_does_not_break_the_json(self):
        """Nor is it cut short: under a UTF-8 locale, macOS `tr` stops at the first byte that
        is not UTF-8 ("Illegal byte sequence"), so every `tr` on HIVE's answer runs in C."""
        self.plugin(self.laptop, "2.9.0")
        for loc in ("C.UTF-8", "en_US.UTF-8"):
            with self.subTest(locale=loc):
                r = subprocess.run(["bash", os.path.join(self.laptop, "scripts",
                                                         "skill_version.sh"), "--check-hive"],
                                   capture_output=True, timeout=60,
                                   env=self.env(FAKE_HIVE_DOWN="latin1", LC_ALL=loc, LANG=loc))
                self.assertEqual(r.returncode, 5, r.stderr)
                self.assertNotIn(b"Illegal byte sequence", r.stderr)
                j = json.loads(r.stdout.decode("utf-8"))  # strict: raises on a stray byte
                self.assertIn("Verbindung zu berprft  fehlgeschlagen", j["say"])

    def test_no_hive_login_skips_silently_without_a_connection(self):
        """Everyone outside the Core, with no HIVE: unaffected."""
        self.plugin(self.laptop, "2.9.0")
        for args in ((), ("--mode", "local")):
            with self.subTest(args=args):
                r = self.run_sv("--check-hive", *args, login=False)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertEqual(r.stderr, "")
                j = json.loads(r.stdout)
                self.assertEqual((j["status"], j["say"], j["next"]), ("skipped", None, None))
        self.assertEqual(self.n_calls(), 0)

    def test_in_hive_remote_a_missing_login_is_said_not_skipped(self):
        """Review: nothing wrote hive.env, so the check skipped itself for the very staff it was
        for. In hive_remote a login is expected; without one, exit 5 and say how to save it."""
        self.plugin(self.laptop, "2.9.0")
        r = self.run_sv("--check-hive", "--mode", "hive_remote", login=False)
        self.assertEqual(r.returncode, 5, r.stderr)
        j = json.loads(r.stdout)
        self.assertEqual(j["status"], "no_login")
        self.assertIn("check_access.sh", j["say"])
        self.assertEqual(self.n_calls(), 0)

    def test_a_login_saved_in_hive_env_is_used(self):
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.9.0")
        cfg = os.path.join(self.d, "hive.env")
        with open(cfg, "w") as fh:
            fh.write(f"HIVE_USER='gabrig'\nHIVE_KEY='{self.key}'\n")
        j = self.check(0, "--mode", "hive_remote", HIVE_ENV_FILE=cfg, HIVE_EXEC=self.fake,
                       login=False)
        self.assertEqual(j["status"], "in_step")

    def test_this_copy_without_plugin_json_cannot_be_compared(self):
        self.plugin(self.laptop, None)
        self.hive_copy("2.9.0")
        j = self.check(5)
        self.assertEqual((j["status"], j["local"]), ("local_unknown", None))
        self.assertEqual(self.n_calls(), 0)

    def test_arguments(self):
        self.plugin(self.laptop, "2.9.0")
        for args in (("extra",), ("--mode",), ("--mode", "cluster")):
            with self.subTest(args=args):
                self.assertEqual(self.run_sv("--check-hive", *args).returncode, 2)
        self.assertEqual(self.n_calls(), 0)

    def test_the_hive_side_read_needs_nothing_from_the_hive_copy(self):
        """A 2.6.0 copy has no skill_version.sh; the reader goes up with the command."""
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.6.0")
        os.remove(os.path.join(self.hive_root, "scripts", "skill_version.sh"))
        self.assertEqual(self.check(3)["hive"], "2.6.0")


class ReleaseTrust(Harness):
    """/quobyte/proteomics-grp is group-writable and not sticky: any member can rename or
    replace skill_release/. CURRENT_VERSION is trusted only while it AND its folder belong to a
    Core admin (scripts/core_admins.txt); otherwise it reads as if nothing were published."""

    def setUp(self):
        super().setUp()
        self.plugin(self.laptop, "2.6.0")
        self.hive_copy("2.6.0")
        self.release("2.9.0\n")

    def test_an_admins_release_is_trusted(self):
        j = self.check(0)
        self.assertEqual((j["release"], j["release_state"], j["behind_release"]),
                         ("2.9.0", "published", True))

    def test_a_members_file_or_folder_is_not(self):
        for who in ({"FAKE_FILE_OWNER": "gabrig"}, {"FAKE_DIR_OWNER": "gabrig"},
                    {"FAKE_FILE_OWNER": "", "FAKE_DIR_OWNER": ""}):   # owner unknown: not trusted
            with self.subTest(**who):
                j = self.check(0, **who)
                self.assertEqual((j["release"], j["release_state"], j["behind_release"], j["say"]),
                                 (None, "untrusted", None, None))

    def test_the_real_stat_reads_the_owner(self):
        """No stand-in: GNU `stat -c %U` on HIVE, BSD `stat -f %Su` here on a Mac."""
        j = self.check(0, SKILL_OWNER_CMD="")
        self.assertEqual(j["release_state"], "published")      # the runner is an admin here
        self.admins(text="# nobody who runs this suite\nbrettsp-not-me\n")
        self.hive_copy("2.6.0")
        self.assertEqual(self.check(0, SKILL_OWNER_CMD="")["release_state"], "untrusted")


class CoreAdmins(unittest.TestCase):
    def admins(self, path=None):
        r = subprocess.run(["bash", "-c", '. "$1"; skill_core_admins "$2"', "bash",
                            os.path.join(SCRIPTS, "skill_version.sh"), path or ""],
                           capture_output=True, text=True, timeout=30,
                           env=job_env(tempfile.gettempdir()))
        return r.stdout.splitlines()

    def test_the_shipped_list(self):
        self.assertEqual(self.admins(), ["brettsp"])

    def test_comments_blank_lines_spaces_and_crlf(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "core_admins.txt")
            with open(p, "w", newline="") as fh:
                fh.write("# admins\r\n\r\n  brettsp  \r\n\t# indented comment\n"
                         "gabrig # staff, 2026\n\n#msalemi\n")
            self.assertEqual(self.admins(p), ["brettsp", "gabrig"])
            self.assertEqual(self.admins(os.path.join(tmp, "missing.txt")), [])

    def test_membership_is_the_whole_name(self):
        r = subprocess.run(["bash", "-c", '. "$1"; for u in brettsp brett rettsp ""; do '
                            'skill_is_core_admin "$u" && echo "$u"; done', "bash",
                            os.path.join(SCRIPTS, "skill_version.sh")],
                           capture_output=True, text=True, timeout=30,
                           env=job_env(tempfile.gettempdir()))
        self.assertEqual(r.stdout.split(), ["brettsp"])


class VersionCompare(unittest.TestCase):
    def cmp(self, a, b):
        with tempfile.TemporaryDirectory() as tmp:
            r = subprocess.run(["bash", "-c", '. "$1"; skill_version_cmp "$2" "$3"', "bash",
                                os.path.join(SCRIPTS, "skill_version.sh"), a, b],
                               capture_output=True, text=True, timeout=30, env=job_env(tmp))
        return r.stdout.strip() if r.returncode == 0 else None

    def test_numbers_of_any_length(self):
        self.assertEqual(self.cmp("2.10.0", "2.9.0"), "1")
        self.assertEqual(self.cmp("2.9", "2.9.0"), "0")
        self.assertEqual(self.cmp("007.1", "7.1"), "0")
        self.assertEqual(self.cmp("1.99999999999999999999999", "1.99999999999999999999998"), "1")
        self.assertEqual(self.cmp("99999999999999999999", "100000000000000000000"), "-1")
        for bad in ("2.9.0-dev", "", "1..2", "v2"):
            self.assertIsNone(self.cmp(bad, "1.0"), bad)


class ThroughTheRealHiveExec(Harness):
    """The multi-line command survives hive_exec.sh's `bash -l -c <quoted>` over ssh."""

    def test_end_to_end(self):
        self.plugin(self.laptop, "2.9.0")
        self.hive_copy("2.6.0")
        self.release("2.9.0\n")
        fakebin = os.path.join(self.d, "bin")
        os.makedirs(fakebin)
        # the remote command is ssh's last argument; run it as HIVE's sshd would
        sh_exe(os.path.join(fakebin, "ssh"),
               'echo call >> "$FAKE_CALLS"\n'
               'for last in "$@"; do :; done\nHOME="$FAKE_HIVE_HOME" exec bash -c "$last"\n')
        env = self.env(PATH=fakebin + os.pathsep + os.environ.get("PATH", ""),
                       HIVE_SSH_MUX="0")
        env.pop("HIVE_EXEC")                                     # the real one, beside it
        r = subprocess.run(["bash", os.path.join(self.laptop, "scripts", "skill_version.sh"),
                            "--check-hive", "--mode", "hive_remote"], capture_output=True,
                           text=True, env=env, timeout=60)
        self.assertEqual(r.returncode, 3, r.stdout + r.stderr)
        j = json.loads(r.stdout)
        self.assertEqual((j["status"], j["hive"], j["release"], j["jobs"]),
                         ("hive_behind", "2.6.0", "2.9.0", 0))
        self.assertEqual(self.n_calls(), 1)


class PublishRelease(Harness):
    def published(self):
        with open(os.path.join(self.release_dir, "CURRENT_VERSION")) as fh:
            return fh.read()

    def test_publish_then_the_check_reads_it(self):
        self.plugin(self.laptop, "2.9.0")
        r = self.run_sv("--publish-release")
        self.assertEqual(r.returncode, 0, r.stderr)
        j = json.loads(r.stdout)
        self.assertEqual((j["published"], j["previous"], j["main_checked"]), ("2.9.0", None, False))
        text = self.published()
        self.assertEqual(text.splitlines()[0], "2.9.0")
        self.assertTrue(all(l.startswith("#") for l in text.splitlines()[1:]), text)
        mode = stat.S_IMODE(os.stat(os.path.join(self.release_dir, "CURRENT_VERSION")).st_mode)
        self.assertTrue(mode & stat.S_IRGRP, oct(mode))
        self.assertEqual([n for n in os.listdir(self.release_dir) if n.startswith(".")], [],
                         "no temporary file left behind")
        # a laptop behind it now hears about it
        self.plugin(self.laptop, "2.8.0")
        self.hive_copy("2.8.0")
        self.assertIs(self.check(0)["behind_release"], True)

    def test_going_back_needs_allow_older(self):
        self.release("2.9.0\n")
        self.plugin(self.laptop, "2.8.0")
        r = self.run_sv("--publish-release")
        self.assertEqual(r.returncode, 2)
        self.assertIn("--allow-older", r.stderr)
        self.assertEqual(self.published(), "2.9.0\n")
        r = self.run_sv("--publish-release", "--allow-older")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(json.loads(r.stdout)["previous"], "2.9.0")
        self.assertEqual(self.published().splitlines()[0], "2.8.0")

    def test_only_a_release_number_is_published(self):
        for v in ("2.9.0-dev", None):
            with self.subTest(v=v):
                self.plugin(self.laptop, v)
                r = self.run_sv("--publish-release")
                self.assertEqual(r.returncode, 2)
                self.assertFalse(os.path.exists(self.release_dir))
        self.assertEqual(self.n_calls(), 0)

    def test_a_v_prefix_is_published_as_the_number(self):
        self.plugin(self.laptop, "v2.9.0")
        r = self.run_sv("--publish-release")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.published().splitlines()[0], "2.9.0")

    def test_it_needs_a_hive_login(self):
        self.plugin(self.laptop, "2.9.0")
        r = self.run_sv("--publish-release", login=False)
        self.assertEqual(r.returncode, 2)
        self.assertIn("HIVE login", r.stderr)

    def test_only_an_admin_publishes(self):
        self.plugin(self.laptop, "2.9.0")
        self.admins()                                  # the shipped list: not the runner
        r = self.run_sv("--publish-release")
        self.assertEqual(r.returncode, 2)
        self.assertIn(f"'{self.me}' is not a Proteomics Core admin", r.stderr)
        self.assertFalse(os.path.exists(self.release_dir))

    def test_a_folder_a_member_made_is_not_published_into(self):
        self.plugin(self.laptop, "2.9.0")
        os.makedirs(self.release_dir)
        r = self.run_sv("--publish-release", FAKE_DIR_OWNER="gabrig")
        self.assertEqual(r.returncode, 2)
        self.assertIn("owned by 'gabrig', not a Core admin", r.stderr)
        self.assertEqual(os.listdir(self.release_dir), [])

    def test_an_account_outside_the_core_cannot_publish(self):
        self.plugin(self.laptop, "2.9.0")
        os.chmod(self.grp, 0o555)
        r = self.run_sv("--publish-release")
        self.assertEqual(r.returncode, 2)
        self.assertIn("Proteomics Core account", r.stderr)

    @unittest.skipUnless(shutil.which("git"), "needs git")
    def test_a_checkout_publishes_only_what_origin_main_ships(self):
        """From a git checkout the version must be the one origin/main has (as last fetched --
        no network); the installed plugin, not a checkout, is not checked (main_checked)."""
        g = ["git", "-C", self.laptop, "-c", "user.name=t", "-c", "user.email=t@t"]
        self.plugin(self.laptop, "2.9.0")
        for cmd in (["init", "-q"], ["add", ".claude-plugin/plugin.json"],
                    ["commit", "-qm", "2.9.0"], ["update-ref", "refs/remotes/origin/main", "HEAD"]):
            subprocess.run(g + cmd, check=True, capture_output=True)
        self.plugin(self.laptop, "2.10.0")                     # bumped, not merged yet
        r = self.run_sv("--publish-release")
        self.assertEqual(r.returncode, 2)
        self.assertIn("origin/main", r.stderr)
        self.assertFalse(os.path.exists(self.release_dir))
        self.plugin(self.laptop, "2.9.0")
        r = self.run_sv("--publish-release")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIs(json.loads(r.stdout)["main_checked"], True)


class Documented(unittest.TestCase):
    def read(self, *rel):
        with open(os.path.join(SKILL, *rel), encoding="utf-8") as fh:
            return fh.read()

    def test_skill_md_step_0_runs_the_check_and_says_what_each_exit_means(self):
        md = self.read("SKILL.md")
        step0 = md[md.index("### 0. One-time setup"):md.index("### 0b.")]
        self.assertIn("bash scripts/skill_version.sh --check-hive --mode hive_remote", step0)
        for code in ("exit 0", "exit 3", "exit 4", "exit 5", "exit 6"):
            self.assertIn(code, step0)
        self.assertIn("--put-skill", step0)

    def test_step_0_puts_the_update_question_first(self):
        """Brett: staff must hear it plainly -- relayed before any other work, an offer to
        update now (the steps that work), the decline recorded, and auto-update recommended."""
        md = self.read("SKILL.md")
        step0 = " ".join(md[md.index("### 0. One-time setup"):md.index("### 0b.")].split())
        for text in ("Relay it FIRST, before any other work",
                     "Your copy of the skill is older than the Core's current release",
                     "ask whether to update now",
                     "claude plugin update ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core",
                     "/reload-plugins",
                     "record that they declined",
                     "always recommend turning on auto-update",
                     "Marketplaces → ucdavis-proteomics-core → Enable auto-update"):
            self.assertIn(text, step0)
        self.assertLess(step0.index("before any other work"), step0.index("**exit 0**"))

    def test_who_publishes_the_release_is_written_down(self):
        access = self.read("references", "access.md")
        self.assertIn("skill_version.sh --publish-release", access)
        self.assertIn("/quobyte/proteomics-grp/skill_release/CURRENT_VERSION", access)
        self.assertIn("scripts/core_admins.txt", access)
        self.assertIn("untrusted", access)

    def test_no_doc_says_marketplace_update_updates_the_plugin(self):
        """`/plugin marketplace update` refreshes the catalogue only (Claude Code 2.1.285)."""
        repo = os.path.dirname(os.path.dirname(SKILL))
        texts = [self.read("references", "access.md"), self.read("SKILL.md")]
        page = os.path.join(repo, "docs", "skill-install.html")
        if os.path.isfile(page):                     # the skill installed on its own has no docs/
            with open(page, encoding="utf-8") as fh:
                texts.append(fh.read())
            self.assertIn("Enable auto-update", texts[-1])
            self.assertIn(UPDATE, texts[-1])
            self.assertNotIn("Updates arrive automatically", texts[-1])
        for t in texts:
            self.assertNotIn('data-copy="/plugin marketplace update', t)
            self.assertNotIn("Claude Code `/plugin marketplace update ucdavis-proteomics-core`,", t)

    # Every place that shows how to install. Brett pasted the README's Install block -- two
    # /plugin commands in one code block -- into Claude Code, and it failed: they have to be
    # entered one at a time. The terminal one-liner (checked on claude 2.1.285) is the easy path.
    INSTALL_DOCS = ("README_GITHUB.md", "README.md", "docs/skill-install.html",
                    "docs/STUDENT_SETUP.md",
                    "skill/ucdavis-proteomics-core-pipeline/SKILL.md",
                    "skill/ucdavis-proteomics-core-pipeline/references/install.md",
                    "skill/ucdavis-proteomics-core-pipeline/references/access.md")
    ONE_LINE = ("claude plugin marketplace add bsphinney/DE-LIMP && "
                "claude plugin install ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core")

    def install_docs(self):
        repo = os.path.dirname(os.path.dirname(SKILL))
        found = {}
        for rel in self.INSTALL_DOCS:
            p = os.path.join(repo, rel)
            if os.path.isfile(p):              # the skill installed on its own has no repo docs
                with open(p, encoding="utf-8") as fh:
                    found[rel] = fh.read()
        return found

    @staticmethod
    def copy_blocks(rel, text):
        """What a reader copies in one go: a fenced ``` block of a Markdown file; a page's
        <pre> block and each Copy button's data-copy."""
        import html
        import re
        if rel.endswith(".html"):
            blocks = re.findall(r"<pre[^>]*>(.*?)</pre>", text, re.S)
            blocks += re.findall(r'data-copy="([^"]*)"', text)
            return [html.unescape(re.sub(r"<[^>]+>", "", b)) for b in blocks]
        return re.findall(r"^[ \t]*```[^\n]*\n(.*?)^[ \t]*```", text, re.S | re.M)

    def test_no_copy_block_holds_two_plugin_commands(self):
        docs = self.install_docs()
        if "README_GITHUB.md" not in docs:
            self.skipTest("no repo docs beside this skill")
        bad = []
        for rel, text in docs.items():
            for b in self.copy_blocks(rel, text):
                cmds = [l for l in b.splitlines() if l.strip().startswith("/plugin")]
                if len(cmds) > 1:
                    bad.append(f"{rel}: {cmds}")
        self.assertEqual(bad, [], "a slash command per block: pasted together they fail")
        # the check sees blocks at all (a regex that matched nothing would pass anything)
        self.assertTrue(any("/plugin install" in b for b in
                            self.copy_blocks("README_GITHUB.md", docs["README_GITHUB.md"])))
        self.assertTrue(any("/plugin install" in b for b in
                            self.copy_blocks("x.html", docs["docs/skill-install.html"])))

    def test_the_one_line_terminal_install_is_offered(self):
        docs = self.install_docs()
        if "README_GITHUB.md" not in docs:
            self.skipTest("no repo docs beside this skill")
        for rel in ("README_GITHUB.md", "docs/skill-install.html", "docs/STUDENT_SETUP.md",
                    "skill/ucdavis-proteomics-core-pipeline/references/install.md"):
            with self.subTest(rel=rel):
                self.assertTrue(any(b.strip() == self.ONE_LINE
                                    for b in self.copy_blocks(rel, docs[rel])), rel)
        self.assertIn("press Enter, then the second", docs["README_GITHUB.md"])
        self.assertIn("Enable\nauto-update", docs["README_GITHUB.md"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
