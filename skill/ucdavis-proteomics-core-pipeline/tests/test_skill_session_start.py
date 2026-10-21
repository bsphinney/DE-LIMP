#!/usr/bin/env python3
"""
One call at the start of a conversation (skill_version.sh --session-start, SKILL.md step 0).

Step 0 had the agent run --check-update, then --check-hive in hive_remote, and follow two trees
of exit codes (0/3/5 and 0/3/4/5/6, with `behind_release` able to come with four of them). The
one call runs both and says what to do in fields: `up_to_date` (the only case in which the skill
may be called up to date), `relay_first`, `ask_update`, `ask`, `next`, `say`; exit 2 when the
user must answer or act first, 3 when `next` is to be run, else 0.

No test uses the network: GitHub main is a plugin.json read through a file:// URL
(SKILL_MAIN_URL), the installed plugin is installed_plugins.json beside a stand-in catalogue
(SKILL_MARKETPLACE_DIR, no `claude` on SKILL_CLAUDE), and HIVE is test_skill_drift's stand-ins
(HIVE_EXEC, SKILL_SQUEUE, SKILL_OWNER_CMD).
"""
import json
import os
import shutil
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from test_skill_drift import Harness as DriftHarness, PUT     # noqa: E402

UPDATE_CMD = "bash scripts/skill_version.sh --update"
ACTIONS = ("up_to_date", "relay_first", "ask_update", "ask", "next", "say")


def older(local, core):
    return (f"Your copy of the skill is older than the Core's current release (you have {local}, "
            f"the Core is on {core}). Updating takes a minute.")


@unittest.skipUnless(shutil.which("curl"), "needs curl (main is read through a file:// URL)")
class Harness(DriftHarness):
    def setUp(self):
        super().setUp()
        self.main_file = os.path.join(self.d, "github_main", "plugin.json")
        self.plugins = os.path.join(self.d, "claude_cfg", "plugins")
        self.market = os.path.join(self.plugins, "marketplaces", "ucdavis-proteomics-core")
        os.makedirs(self.market)

    def main(self, version):
        """GitHub main ships `version`."""
        os.makedirs(os.path.dirname(self.main_file), exist_ok=True)
        with open(self.main_file, "w") as fh:
            json.dump({"name": "ucdavis-proteomics-core-pipeline", "version": version}, fh)

    def installed(self, version):
        """installed_plugins.json beside the catalogue, as Claude Code 2.1.295 writes it."""
        with open(os.path.join(self.plugins, "installed_plugins.json"), "w") as fh:
            json.dump({"version": 2, "plugins": {
                "ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core":
                    [{"scope": "user", "version": version}]}}, fh)

    def given(self, laptop, main, hive=None, release=None):
        self.plugin(self.laptop, laptop)
        self.main(main)
        self.installed(laptop)
        if hive:
            self.hive_copy(hive)
        if release:
            self.release(release + "\n")

    def start(self, want_rc, *args, **extra):
        env = dict(SKILL_MAIN_URL="file://" + self.main_file, SKILL_MARKETPLACE_DIR=self.market,
                   SKILL_CLAUDE=os.path.join(self.d, "no_claude_here"))
        env.update(extra)
        r = self.run_sv("--session-start", *args, **env)
        self.assertEqual(r.returncode, want_rc, r.stdout + r.stderr)
        j = json.loads(r.stdout)
        self.assertEqual(list(j)[:len(ACTIONS)], list(ACTIONS))
        return j

    def remote(self, want_rc, **extra):
        return self.start(want_rc, "--mode", "hive_remote", **extra)


class SessionStart(Harness):
    def test_current_everywhere_says_nothing(self):
        self.given("2.13.1", "2.13.1", hive="2.13.1", release="2.13.1")
        j = self.remote(0)
        self.assertEqual([j[k] for k in ACTIONS], [True, None, False, None, None, None])
        self.assertEqual((j["update"]["status"], j["update"]["exit"]), ("current", 0))
        self.assertEqual((j["hive"]["status"], j["hive"]["exit"]), ("in_step", 0))
        self.assertEqual(self.n_calls(), 1)

    def test_behind_main_relays_first_and_asks_to_update(self):
        self.given("2.12.0", "2.13.1")
        j = self.start(2)
        self.assertEqual(j["relay_first"], older("2.12.0", "2.13.1"))
        self.assertEqual((j["up_to_date"], j["ask_update"], j["ask"], j["next"], j["say"],
                          j["hive"]), (False, True, None, None, None, None))
        self.assertEqual((j["update"]["status"], j["update"]["next"], j["update"]["exit"]),
                         ("behind", UPDATE_CMD, 3))

    def test_the_question_comes_before_putting_the_hive_copy(self):
        """Behind main and HIVE behind this computer: exit 2 (ask first), `next` kept for a
        declined update."""
        self.given("2.12.0", "2.13.1", hive="2.11.2")
        j = self.remote(2)
        self.assertEqual((j["ask_update"], j["next"], j["hive"]["exit"]), (True, PUT, 3))
        self.assertIn("Putting this computer's 2.12.0 there", j["say"])

    def test_reload_needed_waits_for_the_reload(self):
        self.given("2.12.0", "2.13.1", hive="2.11.2")
        self.installed("2.13.1")                    # installed, but this session runs 2.12.0
        j = self.remote(2)
        self.assertIn("/reload-plugins", j["relay_first"])
        self.assertIn("start again", j["relay_first"])
        self.assertEqual((j["up_to_date"], j["ask_update"], j["ask"], j["next"], j["say"],
                          j["hive"]), (False, False, None, None, None, None))
        self.assertEqual((j["update"]["status"], j["update"]["exit"]), ("reload_needed", 3))
        self.assertEqual(self.n_calls(), 0, "HIVE is checked after the reload, not before")

    def test_main_unreadable_is_said_and_never_up_to_date(self):
        self.given("2.13.1", "2.13.1")
        j = self.start(0, SKILL_MAIN_URL="file://" + os.path.join(self.d, "nowhere.json"))
        self.assertEqual((j["up_to_date"], j["relay_first"], j["ask_update"], j["next"]),
                         (False, None, False, None))
        self.assertIn("not a confirmation", j["say"])
        self.assertEqual(j["say"], j["update"]["say"])
        self.assertEqual((j["update"]["status"], j["update"]["exit"]), ("main_unknown", 5))


    def test_an_older_hive_copy_is_put_unasked(self):
        self.given("2.13.1", "2.13.1", hive="2.12.0", release="2.13.1")
        j = self.remote(3)
        self.assertEqual((j["up_to_date"], j["ask_update"], j["ask"], j["next"]),
                         (True, False, None, PUT))
        self.assertIn("HIVE runs skill 2.12.0", j["say"])
        self.assertEqual((j["hive"]["status"], j["hive"]["exit"]), ("hive_behind", 3))

    def test_a_newer_unreleased_hive_build_is_a_question(self):
        self.given("2.13.1", "2.13.1", hive="2.14.0", release="2.13.1")
        j = self.remote(2)
        self.assertEqual((j["up_to_date"], j["relay_first"], j["ask_update"], j["next"],
                          j["say"]), (True, None, False, None, None))
        self.assertIn("newer skill (2.14.0", j["ask"])
        self.assertIn("not released yet", j["ask"])
        self.assertIn(PUT, j["ask"])
        self.assertTrue(j["ask"].endswith("?"), j["ask"])
        self.assertEqual((j["hive"]["status"], j["hive"]["exit"]), ("hive_ahead", 4))
        # no release readable: still a question, never a put
        j = self.remote(2, SKILL_RELEASE_DIR=os.path.join(self.d, "nowhere", "skill_release"))
        self.assertIn("Keep HIVE's 2.14.0", j["ask"])
        self.assertNotIn("not released yet", j["ask"])
        self.assertIsNone(j["next"])

    def test_a_newer_released_hive_copy_is_the_update_question(self):
        """GitHub unreadable from here, HIVE has the release: relay it, ask to update, no
        which-copy question."""
        self.given("2.12.0", "2.13.1", hive="2.13.1", release="2.13.1")
        j = self.remote(2, SKILL_MAIN_URL="file://" + os.path.join(self.d, "nowhere.json"))
        self.assertEqual(j["relay_first"], older("2.12.0", "2.13.1"))
        self.assertEqual((j["up_to_date"], j["ask_update"], j["ask"], j["next"]),
                         (False, True, None, None))
        self.assertIn("not a confirmation", j["say"])
        self.assertNotIn("newer than this computer", j["say"])
        self.assertEqual((j["hive"]["exit"], j["hive"]["behind_release"]), (4, True))

    def test_a_hive_copy_in_use_by_jobs_is_a_question(self):
        self.given("2.13.1", "2.13.1", hive="2.12.0")
        lst = self.jobs(f"python3 {self.hive_root}/scripts/notify_slack.py job-end",
                        "diann --f a.d --out x")
        j = self.remote(2, FAKE_SQUEUE_LIST=lst)
        self.assertEqual((j["up_to_date"], j["ask_update"], j["next"], j["say"]),
                         (True, False, None, None))
        self.assertIn("1 of your 2 SLURM job(s)", j["ask"])
        self.assertIn(f"Replace it now ({PUT}), or after those jobs finish?", j["ask"])
        self.assertEqual((j["hive"]["status"], j["hive"]["exit"]), ("hive_behind", 6))
        j = self.remote(2, FAKE_SQUEUE_FAIL=1)
        self.assertIn("squeue on HIVE answered with an error", j["ask"])
        self.assertIn("once squeue answers", j["ask"])

    def test_behind_the_release_with_hive_in_step(self):
        """main's plugin.json says 2.12.0 (current), but the Core's release is 2.13.1: never
        'up to date', and the release sentence is relay_first, not repeated in `say`."""
        self.given("2.12.0", "2.12.0", hive="2.12.0", release="2.13.1")
        j = self.remote(2)
        self.assertEqual(j["relay_first"], older("2.12.0", "2.13.1"))
        self.assertEqual((j["up_to_date"], j["ask_update"], j["ask"], j["next"], j["say"]),
                         (False, True, None, None, None))
        self.assertEqual((j["update"]["exit"], j["hive"]["exit"], j["hive"]["behind_release"]),
                         (0, 0, True))
        self.assertIn("current release is 2.13.1", j["hive"]["say"])

    def test_behind_both_names_main(self):
        self.given("2.12.0", "2.13.1", hive="2.12.0", release="2.13.0")
        j = self.remote(2)
        self.assertEqual(j["relay_first"], older("2.12.0", "2.13.1"))
        self.assertIsNone(j["say"])

    def test_hives_setup_json_note_is_said(self):
        self.given("2.13.1", "2.13.1", hive="2.12.0", release="2.13.1")
        sj = os.path.join(self.hive_home, ".proteomics-pipeline", "setup.json")
        os.makedirs(os.path.dirname(sj))
        with open(sj, "w") as fh:
            fh.write('{\n  "ready_for": {"de": true},\n  "notes": []\n}\n')
        j = self.remote(3)
        self.assertIn("Putting this computer's 2.13.1 there.", j["say"])
        self.assertIn("setup.sh --check", j["say"])
        # with jobs in the way the status part is the question; the note stays in `say`
        lst = self.jobs(f"python3 {self.hive_root}/scripts/record_run.py")
        j = self.remote(2, FAKE_SQUEUE_LIST=lst)
        self.assertTrue(j["say"].startswith("HIVE's setup.json"), j["say"])
        self.assertIsNotNone(j["ask"])

    def test_hive_not_checked_is_said_in_one_line(self):
        self.given("2.13.1", "2.13.1")
        j = self.remote(0, FAKE_HIVE_DOWN=1)
        self.assertEqual((j["up_to_date"], j["next"], j["ask"], j["hive"]["exit"]),
                         (True, None, None, 5))
        self.assertIn("Could not reach HIVE", j["say"])
        # its words pass through `say` whole, quotes and backslashes included
        down = os.path.join(self.d, "hive_down.sh")
        with open(down, "w") as fh:
            fh.write('#!/usr/bin/env bash\necho \'ssh: key "id_x" in C:\\keys refused\' >&2\n'
                     'exit 255\n')
        os.chmod(down, 0o755)
        j = self.remote(0, HIVE_EXEC=down)
        self.assertIn('(ssh: key "id_x" in C:\\keys refused)', j["say"])
        self.assertEqual(j["say"], j["hive"]["say"])
        j = self.remote(0, login=False)
        self.assertIn("check_access.sh", j["say"])

    def test_without_hive_remote_hive_is_not_checked(self):
        self.given("2.13.1", "2.13.1", hive="2.12.0")
        for args in ((), ("--mode", "local")):
            with self.subTest(args=args):
                j = self.start(0, *args)
                self.assertEqual([j[k] for k in ACTIONS], [True, None, False, None, None, None])
                self.assertIsNone(j["hive"])
        self.assertEqual(self.n_calls(), 0)

    def test_arguments(self):
        self.given("2.13.1", "2.13.1")
        for args in (("x",), ("--mode",), ("--mode", "cluster"), ("--mode", "local", "x")):
            with self.subTest(args=args):
                r = self.run_sv("--session-start", *args)
                self.assertEqual((r.returncode, r.stdout), (2, ""))
                self.assertIn("usage: --session-start", r.stderr)
        self.assertEqual(self.n_calls(), 0)

    def test_help_prints_the_whole_header(self):
        r = self.run_sv("--help")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("--session-start [--mode hive_remote|local]", r.stdout)
        self.assertIn("no test uses the network.", r.stdout)    # the BEHIND GITHUB MAIN part
        self.assertTrue(r.stdout.rstrip().endswith("else 0."), r.stdout[-200:])
        self.assertNotIn("=====", r.stdout)
        self.assertNotIn("SKILL_VERSION_UNKNOWN=", r.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
