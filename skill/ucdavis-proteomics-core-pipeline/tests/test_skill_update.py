#!/usr/bin/env python3
"""
The documented skill update must reach GitHub main on a Windows laptop with no GitHub SSH key,
and must never say "already latest" while the computer is behind main
(skill_version.sh --check-update / --update, SKILL.md step 0).

On 2026-10-07 (msalemi, PROT_0803, Git Bash) `claude plugin update
ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core` warned "marketplace not refreshed: SSH
host key is not in your known_hosts", then "SSH authentication failed", and said "already at the
latest version (2.10.0)" while main shipped 2.11.2: it compared against its own stale catalogue.
The catalogue clone had an HTTPS remote; `git -C <clone> pull --ff-only` and the same update
took it to 2.11.2.

No test uses the network. GitHub is a local bare repository, main's plugin.json is read through
a file:// URL (SKILL_MAIN_URL), `claude` is a stand-in (SKILL_CLAUDE) that, like the one on the
laptop, never refreshes its catalogue itself, and `ssh` on PATH is a stand-in that records any
call and fails, so a test sees if anything tried SSH.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
SV = os.path.join(SCRIPTS, "skill_version.sh")
sys.path.insert(0, HERE)
from job_env import job_env     # noqa: E402

ID = "ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core"
UPDATE = "claude plugin update " + ID
UPDATE_CMD = "bash scripts/skill_version.sh --update"
PLUGIN_REL = os.path.join("skill", "ucdavis-proteomics-core-pipeline", ".claude-plugin",
                          "plugin.json")

# The stand-in `claude`: the installed version lives in installed_plugins.json beside the
# catalogue, as Claude Code 2.1.295 keeps it. `plugin update` installs what the catalogue on
# disk offers and never refreshes it; FAKE_CLAUDE_SSH=1 adds the laptop's SSH warnings.
FAKE_CLAUDE = r'''#!/usr/bin/env bash
echo "$* PREFER_HTTPS=${CLAUDE_CODE_PLUGIN_PREFER_HTTPS:-}" >> "$FAKE_CLAUDE_LOG"
ID="ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core"
STATE="$FAKE_PLUGINS/installed_plugins.json"
cat_v() { sed -nE 's/.*"version"[[:space:]]*:[[:space:]]*"([^"]*)".*/\1/p' \
  "$FAKE_PLUGINS/marketplaces/ucdavis-proteomics-core/skill/ucdavis-proteomics-core-pipeline/.claude-plugin/plugin.json" | head -n1; }
inst_v() { [ -f "$STATE" ] && sed -nE 's/.*"version": "([^"]*)".*/\1/p' "$STATE" | tail -n1; }
write_state() {
  [ "${FAKE_CLAUDE_NO_FILE:-0}" = 1 ] && { echo "$1" > "$FAKE_PLUGINS/.hidden_state"; return; }
  printf '{\n  "version": 2,\n  "plugins": {\n    "%s": [\n      {\n        "scope": "user",\n        "installPath": "%s/cache/ucdavis-proteomics-core/ucdavis-proteomics-core-pipeline/%s",\n        "version": "%s"\n      }\n    ]\n  }\n}\n' "$ID" "$FAKE_PLUGINS" "$1" "$1" > "$STATE"
}
case "$1 $2" in
  "plugin list")
    if [ "${FAKE_CLAUDE_NO_JSON:-0}" = 1 ]; then echo "error: unknown option '--json'" >&2; exit 1; fi
    v="$(inst_v)"
    printf '[\n  {\n    "id": "other-plugin@claude-plugins-official",\n    "version": "9.9.9",\n    "scope": "user"\n  }'
    [ -n "$v" ] && printf ',\n  {\n    "id": "%s",\n    "version": "%s",\n    "scope": "user",\n    "enabled": true\n  }' "$ID" "$v"
    printf '\n]\n'; exit 0 ;;
  "plugin update")
    [ "$3" = "$ID" ] || { echo "Plugin \"$3\" not found" >&2; exit 1; }
    echo "Checking for updates for plugin \"$ID\"…"
    if [ "${FAKE_CLAUDE_SSH:-0}" = 1 ]; then
      echo "Warning: marketplace not refreshed: SSH host key is not in your known_hosts" >&2
      echo "SSH authentication failed" >&2
    fi
    [ "${FAKE_CLAUDE_FAIL:-0}" = 1 ] && { echo "Failed to update plugin \"$ID\": boom" >&2; exit 1; }
    c="$(cat_v 2>/dev/null)"; i="$(inst_v)"
    [ -n "$c" ] || { echo "Marketplace \"ucdavis-proteomics-core\" not found" >&2; exit 1; }
    if [ "$c" = "$i" ]; then
      echo "ucdavis-proteomics-core-pipeline is already at the latest version ($i)."; exit 0
    fi
    write_state "$c"
    echo "Plugin \"ucdavis-proteomics-core-pipeline\" updated from $i to $c. Restart to apply."; exit 0 ;;
esac
echo "fake claude: unexpected: $*" >&2; exit 2
'''


def sh_exe(path, body):
    with open(path, "w") as fh:
        fh.write(body if body.startswith("#!") else "#!/usr/bin/env bash\n" + body)
    os.chmod(path, 0o755)
    return path


def plugin_json(version, name="ucdavis-proteomics-core-pipeline"):
    return json.dumps({"name": name, "version": version, "description": "x"}, indent=2) + "\n"


class Harness(unittest.TestCase):
    """A laptop: ~/.claude/plugins with the plugin cache (the copy that runs), the marketplace
    clone of a local "GitHub", and installed_plugins.json."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = os.path.realpath(self._tmp.name)
        self.plugins = os.path.join(self.d, "claude_cfg", "plugins")
        self.market = os.path.join(self.plugins, "marketplaces", "ucdavis-proteomics-core")
        self.claude_log = os.path.join(self.d, "claude.log")
        self.ssh_log = os.path.join(self.d, "ssh.log")
        self.bin = os.path.join(self.d, "bin")
        os.makedirs(self.bin)
        os.makedirs(os.path.join(self.d, "home"))
        sh_exe(os.path.join(self.bin, "ssh"), 'echo "$*" >> "$FAKE_SSH_LOG"\n'
               'echo "Permission denied (publickey)." >&2; exit 255\n')
        self.claude = sh_exe(os.path.join(self.d, "fake_claude"), FAKE_CLAUDE)
        self.github = os.path.join(self.d, "github", "DE-LIMP.git")      # the bare "GitHub"
        self.work = os.path.join(self.d, "maintainer")                    # pushes releases
        self.main_url = "file://" + os.path.join(self.work, PLUGIN_REL)   # raw main plugin.json

    def tearDown(self):
        self._tmp.cleanup()

    # ------------------------------------------------------------------ fixtures
    def git(self, *args, cwd=None):
        r = subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t",
                            "-c", "init.defaultBranch=main", *args], cwd=cwd,
                           capture_output=True, text=True, env=self.git_env(), timeout=60)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return r.stdout.strip()

    def git_env(self):
        e = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
        e.update(HOME=os.path.join(self.d, "home"), GIT_CONFIG_NOSYSTEM="1",
                 PATH=self.bin + os.pathsep + os.environ.get("PATH", ""),
                 FAKE_SSH_LOG=self.ssh_log)
        return e

    def release(self, version):
        """Main ships `version`: committed in the maintainer's checkout and pushed to GitHub."""
        if not os.path.isdir(self.github):
            os.makedirs(self.github)
            self.git("init", "-q", "--bare", self.github)
            self.git("--git-dir", self.github, "symbolic-ref", "HEAD", "refs/heads/main")
            self.git("clone", "-q", self.github, self.work)
            self.git("-C", self.work, "checkout", "-q", "-b", "main")
        p = os.path.join(self.work, PLUGIN_REL)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            fh.write(plugin_json(version))
        self.git("-C", self.work, "add", PLUGIN_REL)
        self.git("-C", self.work, "commit", "-q", "--allow-empty", "-m", version)
        self.git("-C", self.work, "push", "-q", "origin", "HEAD:main")

    def catalogue(self):
        """The marketplace clone, as Claude Code makes it: shallow, of GitHub's main."""
        os.makedirs(os.path.dirname(self.market), exist_ok=True)
        self.git("clone", "-q", "--depth", "1", "--branch", "main", "file://" + self.github,
                 self.market)

    def catalogue_version(self):
        with open(os.path.join(self.market, PLUGIN_REL)) as fh:
            return json.load(fh)["version"]

    def installed(self, version):
        """installed_plugins.json, and the plugin cache copy of `version` (the one that runs)."""
        with open(os.path.join(self.plugins, "installed_plugins.json"), "w") as fh:
            json.dump({"version": 2, "plugins": {ID: [{"scope": "user", "version": version}]}},
                      fh, indent=2)
        return self.running_copy(version)

    def running_copy(self, version):
        root = os.path.join(self.plugins, "cache", "ucdavis-proteomics-core",
                            "ucdavis-proteomics-core-pipeline", version)
        os.makedirs(os.path.join(root, "scripts"), exist_ok=True)
        os.makedirs(os.path.join(root, ".claude-plugin"), exist_ok=True)
        shutil.copy(SV, os.path.join(root, "scripts"))
        with open(os.path.join(root, ".claude-plugin", "plugin.json"), "w") as fh:
            fh.write(plugin_json(version))
        self.copy = root
        return root

    def laptop(self, main="2.11.2", installed="2.10.0"):
        """msalemi's laptop on 2026-10-07: main at 2.11.2, installed and catalogue at 2.10.0."""
        self.release(installed)
        self.catalogue()
        self.release(main)
        self.installed(installed)

    def env(self, **extra):
        e = job_env(self.d, HOME=os.path.join(self.d, "home"), TMPDIR=self.d,
                    PATH=self.bin + os.pathsep + os.environ.get("PATH", ""),
                    SKILL_MAIN_URL=self.main_url, SKILL_CLAUDE=self.claude,
                    FAKE_CLAUDE_LOG=self.claude_log, FAKE_PLUGINS=self.plugins,
                    FAKE_SSH_LOG=self.ssh_log, GIT_CONFIG_NOSYSTEM="1",
                    LC_ALL="C.UTF-8", LANG="C.UTF-8")
        for k in [k for k in e if k.startswith("GIT_SSH") or k in ("CLAUDE_CONFIG_DIR",
                                                                   "SKILL_MARKETPLACE_DIR",
                                                                   "SKILL_REPO_URL")]:
            e.pop(k)
        e.update({k: str(v) for k, v in extra.items()})
        return e

    def run_sv(self, *args, **extra):
        return subprocess.run(["bash", os.path.join(self.copy, "scripts", "skill_version.sh"),
                               *args], capture_output=True, text=True, timeout=120,
                              env=self.env(**extra))

    def sv(self, want_rc, *args, **extra):
        r = self.run_sv(*args, **extra)
        self.assertEqual(r.returncode, want_rc, r.stdout + r.stderr)
        return json.loads(r.stdout)

    def claude_calls(self):
        if not os.path.exists(self.claude_log):
            return []
        with open(self.claude_log) as fh:
            return fh.read().splitlines()

    def updates(self):
        return [c for c in self.claude_calls() if c.startswith("plugin update")]

    def ssh_calls(self):
        if not os.path.exists(self.ssh_log):
            return []
        with open(self.ssh_log) as fh:
            return fh.read().splitlines()

    def assertNeverCurrent(self, say):
        """What staff hear while behind: never 'latest', never 'up to date'."""
        low = (say or "").lower()
        for word in ("latest", "up to date", "up-to-date"):
            self.assertNotIn(word, low, say)


@unittest.skipUnless(shutil.which("git") and shutil.which("curl"), "needs git and curl")
class CheckUpdate(Harness):
    def test_in_step_with_main_says_nothing(self):
        self.laptop(main="2.11.2", installed="2.11.2")
        j = self.sv(0, "--check-update")
        self.assertEqual((j["status"], j["local"], j["main"], j["behind_main"], j["next"],
                          j["say"]), ("current", "2.11.2", "2.11.2", False, None, None))
        self.assertEqual(self.claude_calls(), [], "nothing asked of claude when current")

    def test_a_test_build_ahead_of_main_is_fine(self):
        self.laptop(main="2.11.2", installed="2.11.2")
        self.running_copy("2.12.0-dev")
        self.assertEqual(self.sv(0, "--check-update")["status"], "current")

    def test_behind_main_is_said_with_the_update_that_works(self):
        self.laptop()
        j = self.sv(3, "--check-update")
        self.assertEqual((j["status"], j["local"], j["main"], j["installed"], j["behind_main"],
                          j["next"]), ("behind", "2.10.0", "2.11.2", "2.10.0", True, UPDATE_CMD))
        self.assertIn("2.10.0", j["say"])
        self.assertIn("GitHub main ships 2.11.2", j["say"])
        self.assertIn(UPDATE_CMD, j["say"])
        self.assertIn(UPDATE, j["say"])
        self.assertIn("HTTPS", j["say"])
        self.assertNeverCurrent(j["say"])
        # the by-hand route: the one that worked on the laptop
        self.assertIn("pull --ff-only && " + UPDATE, j["manual"])
        self.assertIn(self.market, j["manual"].replace("\\", ""))
        self.assertEqual(self.updates(), [], "--check-update changes nothing")
        self.assertEqual(self.catalogue_version(), "2.10.0")

    def test_numbers_compare_as_numbers(self):
        self.laptop(main="2.10.0", installed="2.9.0")
        self.assertEqual(self.sv(3, "--check-update")["status"], "behind")    # 2.9 < 2.10
        self.running_copy("2.10.1")
        self.assertEqual(self.sv(0, "--check-update")["status"], "current")

    def test_installed_already_but_this_session_runs_the_old_copy(self):
        self.laptop()
        self.installed("2.11.2")
        self.running_copy("2.10.0")                    # the session still runs 2.10.0
        j = self.sv(3, "--check-update")
        self.assertEqual((j["status"], j["installed"], j["behind_main"], j["next"]),
                         ("reload_needed", "2.11.2", True, None))
        self.assertIn("/reload-plugins", j["say"])

    def test_main_unreadable_is_never_called_up_to_date(self):
        self.laptop()
        j = self.sv(5, "--check-update",
                    SKILL_MAIN_URL="file://" + os.path.join(self.d, "nowhere", "plugin.json"))
        self.assertEqual((j["status"], j["main"], j["behind_main"]), ("main_unknown", None, None))
        self.assertIn("not a confirmation", j["say"])
        self.assertNeverCurrent(j["say"])
        garbled = os.path.join(self.d, "garbled.json")
        with open(garbled, "w") as fh:
            fh.write('{"name": "ucdavis-proteomics-core-pipeline"}\n')
        j = self.sv(5, "--check-update", SKILL_MAIN_URL="file://" + garbled)
        self.assertIn("no readable version", j["say"])

    def test_a_copy_without_plugin_json_cannot_be_compared(self):
        self.laptop()
        os.remove(os.path.join(self.copy, ".claude-plugin", "plugin.json"))
        j = self.sv(5, "--check-update")
        self.assertEqual(j["status"], "local_unknown")

    def test_arguments(self):
        self.laptop()
        for args in (("--check-update", "x"), ("--update", "now")):
            with self.subTest(args=args):
                self.assertEqual(self.run_sv(*args).returncode, 2)
        self.assertEqual(self.updates(), [])


@unittest.skipUnless(shutil.which("git") and shutil.which("curl"), "needs git and curl")
class UpdateOverHttps(Harness):
    def test_the_fake_reproduces_the_laptop(self):
        """Without the refresh, `claude plugin update` alone says 'already at the latest
        version (2.10.0)' with main at 2.11.2 -- what msalemi saw."""
        self.laptop()
        r = subprocess.run([self.claude, "plugin", "update", ID], capture_output=True, text=True,
                           env=self.env(FAKE_CLAUDE_SSH=1), timeout=30)
        self.assertIn("already at the latest version (2.10.0)", r.stdout)
        self.assertIn("marketplace not refreshed", r.stderr)

    def test_the_windows_case_reaches_main(self):
        """Catalogue 2.10.0 behind GitHub's 2.11.2, claude's own refresh failing over SSH: the
        HTTPS pull, then the update, then the version that landed is checked."""
        self.laptop()
        j = self.sv(0, "--update", FAKE_CLAUDE_SSH=1)
        self.assertEqual((j["status"], j["local"], j["main"], j["installed_before"],
                          j["installed"]), ("updated", "2.10.0", "2.11.2", "2.10.0", "2.11.2"))
        self.assertEqual((j["catalogue_before"], j["catalogue"], j["refresh"], j["claude_exit"]),
                         ("2.10.0", "2.11.2", "ok", 0))
        self.assertIn("/reload-plugins", j["say"])
        self.assertEqual(self.catalogue_version(), "2.11.2")
        self.assertEqual(self.updates(), [f"plugin update {ID} PREFER_HTTPS=1"])
        self.assertEqual(self.ssh_calls(), [], "nothing went over SSH")
        # the catalogue is still a fast-forward of GitHub's main, nothing rewritten
        self.assertEqual(self.git("-C", self.market, "rev-parse", "HEAD"),
                         self.git("--git-dir", self.github, "rev-parse", "main"))
        # and now the check agrees, once the new copy runs
        self.running_copy("2.11.2")
        self.assertEqual(self.sv(0, "--check-update")["status"], "current")

    def test_a_failed_refresh_is_reported_as_not_updated_never_latest(self):
        self.laptop()
        self.git("-C", self.market, "remote", "set-url", "origin",
                 "file://" + os.path.join(self.d, "gone", "DE-LIMP.git"))
        j = self.sv(1, "--update", FAKE_CLAUDE_SSH=1)
        self.assertEqual((j["status"], j["installed"], j["refresh"]),
                         ("not_updated", "2.10.0", "failed"))
        self.assertTrue(j["say"].startswith("The update did NOT happen"), j["say"])
        self.assertIn("still has skill 2.10.0", j["say"])
        self.assertIn("GitHub main ships 2.11.2", j["say"])
        self.assertIn("Refreshing the plugin catalogue over HTTPS failed", j["say"])
        self.assertIn("could not refresh its own catalogue", j["say"])
        self.assertIn("pull --ff-only && " + UPDATE, j["say"])
        self.assertNeverCurrent(j["say"])
        self.assertIn("already at the latest version (2.10.0)", j["claude_said"],
                      "claude's own words are kept, for the record")

    def test_claude_saying_latest_on_a_stale_catalogue_is_not_updated(self):
        """The catalogue cannot be refreshed and claude raises no warning: still behind main,
        so still 'did not happen' -- whatever claude printed."""
        self.laptop()
        self.git("-C", self.market, "remote", "set-url", "origin",
                 "file://" + os.path.join(self.d, "gone", "DE-LIMP.git"))
        j = self.sv(1, "--update")
        self.assertEqual(j["status"], "not_updated")
        self.assertIn("compared against its own catalogue, which is at 2.10.0", j["say"])
        self.assertNeverCurrent(j["say"])

    def test_an_ssh_remote_is_pulled_over_https(self):
        """A catalogue whose remote is SSH is pulled from the repository's HTTPS URL instead
        (SKILL_REPO_URL; here the local "GitHub"), never over SSH."""
        self.laptop()
        self.git("-C", self.market, "remote", "set-url", "origin",
                 "git@github.com:bsphinney/DE-LIMP.git")
        url = "file://" + self.github
        j = self.sv(0, "--update", SKILL_REPO_URL=url)
        self.assertEqual((j["status"], j["installed"], j["refresh"]),
                         ("updated", "2.11.2", "ok"))
        self.assertEqual(self.ssh_calls(), [])
        self.assertIn(f"pull --ff-only {url} main && ", j["manual"])

    def test_nothing_is_run_when_already_current(self):
        self.laptop(main="2.11.2", installed="2.11.2")
        j = self.sv(0, "--update")
        self.assertEqual((j["status"], j["refresh"], j["claude_exit"]),
                         ("current", "not_needed", None))
        self.assertEqual(self.claude_calls(), [])

    def test_installed_but_not_reloaded_is_not_updated_again(self):
        self.laptop()
        self.installed("2.11.2")
        self.running_copy("2.10.0")
        j = self.sv(0, "--update")
        self.assertEqual((j["status"], j["installed"]), ("reload_needed", "2.11.2"))
        self.assertEqual(self.updates(), [])
        self.assertIn("/reload-plugins", j["say"])

    def test_a_current_catalogue_is_not_pulled(self):
        """The catalogue already has main's version (someone refreshed it): only the update."""
        self.laptop()
        self.git("-C", self.market, "pull", "-q", "--ff-only")
        head = self.git("-C", self.market, "rev-parse", "HEAD")
        self.git("-C", self.market, "remote", "set-url", "origin",
                 "file://" + os.path.join(self.d, "gone", "DE-LIMP.git"))  # a pull would fail
        j = self.sv(0, "--update")
        self.assertEqual((j["status"], j["refresh"], j["installed"]),
                         ("updated", "not_needed", "2.11.2"))
        self.assertEqual(self.git("-C", self.market, "rev-parse", "HEAD"), head)

    def test_main_unreadable_still_updates_but_is_not_confirmed(self):
        self.laptop()
        j = self.sv(5, "--update",
                    SKILL_MAIN_URL="file://" + os.path.join(self.d, "nowhere.json"))
        self.assertEqual((j["status"], j["main"], j["refresh"], j["installed"]),
                         ("unverified", None, "ok", "2.11.2"))
        self.assertIn("not a confirmation", j["say"])
        self.assertNeverCurrent(j["say"])

    def test_the_installed_file_is_read_when_list_json_is_missing(self):
        """An older claude without `plugin list --json`: installed_plugins.json beside the
        catalogue says what landed."""
        self.laptop()
        j = self.sv(0, "--update", FAKE_CLAUDE_NO_JSON=1)
        self.assertEqual((j["status"], j["installed_before"], j["installed"]),
                         ("updated", "2.10.0", "2.11.2"))

    def test_an_installed_version_that_cannot_be_read_is_not_confirmed(self):
        self.laptop()
        os.remove(os.path.join(self.plugins, "installed_plugins.json"))
        j = self.sv(5, "--update", FAKE_CLAUDE_NO_JSON=1, FAKE_CLAUDE_NO_FILE=1)
        self.assertEqual((j["status"], j["installed"]), ("unverified", None))
        self.assertIn("claude plugin list", j["say"])
        self.assertNeverCurrent(j["say"])

    def test_a_failing_claude_is_not_updated(self):
        self.laptop()
        j = self.sv(1, "--update", FAKE_CLAUDE_FAIL=1)
        self.assertEqual((j["status"], j["claude_exit"], j["refresh"]),
                         ("not_updated", 1, "ok"))
        self.assertIn("claude plugin update failed (exit 1)", j["say"])
        self.assertIn("boom", j["claude_said"])

    def test_no_claude_on_path_is_not_updated(self):
        self.laptop()
        j = self.sv(1, "--update", SKILL_CLAUDE=os.path.join(self.d, "no_such_claude"))
        self.assertEqual((j["status"], j["claude_exit"], j["installed"]),
                         ("not_updated", None, "2.10.0"))
        self.assertIn("not on PATH", j["say"])
        self.assertIn("Update now", j["say"])

    def test_a_catalogue_with_local_changes_is_never_reset(self):
        """Fast-forward only: a hand-edited catalogue fails the pull and keeps its edit."""
        self.laptop()
        p = os.path.join(self.market, PLUGIN_REL)
        with open(p, "a") as fh:
            fh.write("\n")
        with open(p) as fh:
            edited = fh.read()
        j = self.sv(1, "--update")
        self.assertEqual((j["status"], j["refresh"]), ("not_updated", "failed"))
        self.assertTrue(j["refresh_error"], j)
        with open(p) as fh:
            self.assertEqual(fh.read(), edited)

    def test_no_catalogue_at_all(self):
        self.laptop()
        shutil.rmtree(self.market)
        j = self.sv(1, "--update")
        self.assertEqual((j["status"], j["refresh"]), ("not_updated", "no_clone"))
        self.assertIn("no plugin catalogue", j["say"])

    def test_the_catalogue_is_found_beside_the_running_plugin_cache(self):
        """No CLAUDE_CONFIG_DIR and a HOME elsewhere: the catalogue is the one beside
        plugins/cache/ this copy runs from."""
        self.laptop()
        j = self.sv(3, "--check-update")
        self.assertIn(self.market, j["manual"].replace("\\", ""))
        self.assertFalse(self.market.startswith(os.path.join(self.d, "home")))


class InstalledVersionReader(unittest.TestCase):
    """skill_installed_version: the version `claude plugin list --json` or installed_plugins.json
    gives this plugin -- the lowest of several installs, so one still behind is never current."""

    def read(self, text):
        with tempfile.TemporaryDirectory() as tmp:
            r = subprocess.run(["bash", "-c", '. "$1"; skill_installed_version', "bash", SV],
                               input=text, capture_output=True, text=True, timeout=30,
                               env=job_env(tmp))
        self.assertEqual(r.returncode, 0, r.stderr)
        return r.stdout.strip()

    def test_claude_plugin_list_json(self):
        text = json.dumps([{"id": "code-review@claude-plugins-official", "version": "9.9.9"},
                           {"id": ID, "version": "2.11.2", "scope": "user", "enabled": True}],
                          indent=2)
        self.assertEqual(self.read(text), "2.11.2")

    def test_installed_plugins_json(self):
        text = json.dumps({"version": 2, "plugins": {
            "other@x": [{"scope": "user", "version": "8.0.0"}],
            ID: [{"scope": "user", "installPath": "/c/Users/m/.claude/plugins/cache/"
                  "ucdavis-proteomics-core/ucdavis-proteomics-core-pipeline/2.10.0",
                  "version": "2.10.0", "gitCommitSha": "abc"}]}}, indent=2)
        self.assertEqual(self.read(text), "2.10.0")

    def test_one_line_and_crlf(self):
        self.assertEqual(self.read('[{"id":"%s","version":"2.9.0"}]' % ID), "2.9.0")
        self.assertEqual(self.read('[\r\n {\r\n  "id": "%s",\r\n  "version": "2.9.1"\r\n }\r\n]'
                                   % ID), "2.9.1")

    def test_the_lowest_install_wins(self):
        text = json.dumps([{"id": ID, "version": "2.11.2", "scope": "user"},
                           {"id": ID, "scope": "project", "version": "2.10.0"}])
        self.assertEqual(self.read(text), "2.10.0")

    def test_lookalikes_and_non_versions_do_not_count(self):
        text = json.dumps([{"id": "x-" + ID, "version": "1.0.0"},
                           {"id": ID + "-dev", "version": "1.0.0"},
                           {"id": "ucdavis-proteomics-core-pipeline@elsewhere", "version": "1.0.0"},
                           {"id": ID, "version": "315c4e48967d"}])
        self.assertEqual(self.read(text), "")
        self.assertEqual(self.read(""), "")
        self.assertEqual(self.read("[]"), "")


class UpdateState(unittest.TestCase):
    def state(self, *args):
        with tempfile.TemporaryDirectory() as tmp:
            r = subprocess.run(["bash", "-c", '. "$1"; shift; skill_update_state "$@"', "bash",
                                SV, *args], capture_output=True, text=True, timeout=30,
                               env=job_env(tmp))
        self.assertEqual(r.returncode, 0, r.stderr)
        return r.stdout.strip()

    def test_the_table(self):
        cases = [(("2.10.0", "2.11.2"), "behind"),
                 (("2.9.0", "2.10.0"), "behind"),           # numbers, not text
                 (("2.11.2", "2.11.2"), "current"),
                 (("2.12.0", "2.11.2"), "current"),          # a test build ahead of main
                 (("2.11.2-dev", "2.11.2"), "current"),      # a build of main's release
                 (("v2.11.2", "2.11.2"), "current"),
                 (("2.10.0", "2.11.2", "2.11.2"), "reload_needed"),
                 (("2.10.0", "2.11.2", "2.12.0"), "reload_needed"),
                 (("2.10.0", "2.11.2", "2.10.0"), "behind"),
                 (("2.10.0", "2.11.2", ""), "behind"),
                 (("2.10.0", ""), "unknown"),                # main not read: never "current"
                 (("(unknown — plugin.json not found)", "2.11.2"), "unknown"),
                 (("2.10.0", "next week"), "unknown")]
        for args, want in cases:
            with self.subTest(args=args):
                self.assertEqual(self.state(*args), want)


class Documented(unittest.TestCase):
    def read(self, *rel):
        with open(os.path.join(SKILL, *rel), encoding="utf-8") as fh:
            return fh.read()

    def step0(self):
        md = self.read("SKILL.md")
        return " ".join(md[md.index("### 0. One-time setup"):md.index("### 0b.")].split())

    def test_step_0_checks_main_first_in_every_mode(self):
        s = self.step0()
        self.assertIn("bash scripts/skill_version.sh --check-update", s)
        self.assertLess(s.index("--check-update"),
                        s.index("bash scripts/skill_version.sh --check-hive --mode hive_remote"))
        self.assertIn("every mode", s.lower())

    def test_step_0_updates_through_the_https_route_and_checks_it(self):
        s = self.step0()
        for text in ("bash scripts/skill_version.sh --update",
                     "git -C ~/.claude/plugins/marketplaces/ucdavis-proteomics-core pull --ff-only",
                     "marketplace not refreshed",
                     "The update did NOT happen",
                     "already at the latest version",
                     "CLAUDE_CODE_PLUGIN_PREFER_HTTPS=1"):
            self.assertIn(text, s)
        self.assertIn("Never tell the user the skill is up to date", s)

    def test_the_install_docs_describe_the_refresh_route(self):
        for rel in (("references", "install.md"), ("references", "access.md")):
            with self.subTest(rel=rel):
                t = " ".join(self.read(*rel).split())
                self.assertIn("git -C ~/.claude/plugins/marketplaces/ucdavis-proteomics-core "
                              "pull --ff-only", t)
                self.assertIn("https://github.com/bsphinney/DE-LIMP.git", t)
                self.assertIn("skill_version.sh --update", t)


if __name__ == "__main__":
    unittest.main(verbosity=2)
