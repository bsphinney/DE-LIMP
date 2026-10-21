#!/usr/bin/env python3
"""core_sync.sh: the ONE place the Core flow's pull/push lists live (SKILL.md step 1c.7).

The step used to carry two hand-typed loops of hive_exec.sh calls that the agent copied; a
missed or mistyped name was a report missing a section, or a deliverable missing a file. The
lists moved into scripts/core_sync.sh. HIVE here is a directory tree: HIVE_EXEC is a stand-in
that implements hive_exec.sh's --get/--put by copying within it (exiting as rsync does for a
missing source, 23) and runs a plain command there as HIVE's login shell would.
"""
import json
import os
import shutil
import stat
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "scripts", "core_sync.sh")

GET_ITEMS = ["output/tables", "output/figures", "output/search/report.parquet", "output/AUDIT.md",
             "output/SAMPLE_QUALITY.md", "output/methods.md", "session.json",
             "input/submission.json", "input/samples.tsv", "input/raw_files.txt",
             "input/conditions.csv", "input/search.fasta.meta.json", "input/sample_type.json",
             "input/experiment_type.json", "input/locate_decisions.json"]
PUT_ITEMS = ["output/AI_Analysis_Report.md", "output/Analysis_Report.html",
             "output/Analysis_Report.md", "output/Analysis_Report.pdf",
             "output/Analysis_Report_with_audio.html", "output/methods.md", "output/methods.docx",
             "output/podcast"]

# hive_exec.sh stand-in. --get/--put copy within $FAKE_REMOTE (HIVE's / and ~ at once) and, like
# rsync, create only the LAST folder of a destination. FAKE_FAIL=<text>: a transfer naming it
# dies as a dropped connection does (255); FAKE_DOWN=1: every call does.
FAKE_HIVE_EXEC = r'''#!/usr/bin/env bash
printf '%s\n' "$*" >> "$FAKE_LOG"
r() { case "$1" in "~") printf '%s\n' "$FAKE_REMOTE" ;; "~/"*) printf '%s\n' "$FAKE_REMOTE/${1#\~/}" ;; *) printf '%s\n' "$1" ;; esac; }
down() { echo "kex_exchange_identification: read: Operation timed out" >&2; exit 255; }
[ -z "${FAKE_DOWN:-}" ] || down
case "$1" in
  --get|--put)
    [ $# -eq 3 ] || exit 2
    if [ "$1" = --get ]; then src="$(r "$2")"; dst="$3"; else src="$2"; dst="$(r "$3")"; fi
    if [ -n "${FAKE_FAIL:-}" ] && [[ "$2" == *"$FAKE_FAIL"* ]]; then down; fi
    if [ ! -e "$src" ]; then
      echo "rsync: link_stat \"$2\" failed: No such file or directory (2)" >&2
      echo "rsync error: some files/attrs were not transferred (see previous errors) (code 23)" >&2
      exit 23
    fi
    [ -d "$dst" ] || mkdir "${dst%/}" || exit 11
    cp -R "$src" "$dst" ;;
  *)
    echo "Welcome to HIVE (a login shell may print anything first)"
    out="$(cd "$FAKE_REMOTE" && HOME="$FAKE_REMOTE" bash -c "$1")"; rc=$?
    if [ -n "${FAKE_CRLF:-}" ]; then printf '%s\n' "$out" | awk '{ printf "%s\r\n", $0 }'
    else printf '%s\n' "$out"; fi
    exit $rc ;;
esac
'''


def bashes():
    """Every distinct bash here: the one on PATH (bash 5 on CI) and macOS's /bin/bash 3.2."""
    found = []
    for b in (shutil.which("bash"), "/bin/bash"):
        if b and os.path.exists(b) and os.path.realpath(b) not in [os.path.realpath(x) for x in found]:
            found.append(b)
    return found


class CoreSync(unittest.TestCase):
    SESSION = "~/core/PROT_0807/session_1"

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = os.path.realpath(self._tmp.name)
        self.home = self.mk("home")
        self.fake = os.path.join(self.d, "fake_hive_exec.sh")
        with open(self.fake, "w") as fh:
            fh.write(FAKE_HIVE_EXEC)
        os.chmod(self.fake, 0o755)
        self.layout("run")

    def layout(self, name):
        """A fresh HIVE tree, local folder and call log under <tmp>/<name>."""
        self.remote = self.mk(name, "hive")
        self.log = os.path.join(self.d, name, "hive_exec.log")
        self.rsess = os.path.join(self.remote, "core", "PROT_0807", "session_1")
        # a space on purpose: everything is quoted on this side
        self.local = os.path.join(self.d, name, "my core", "session")

    def tearDown(self):
        self._tmp.cleanup()

    def mk(self, *parts):
        p = os.path.join(self.d, *parts)
        os.makedirs(p, exist_ok=True)
        return p

    def write(self, root, rel, text="x"):
        p = os.path.join(root, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def full_session(self):
        for rel in GET_ITEMS:
            if rel in ("output/tables", "output/figures"):
                self.write(self.rsess, rel + "/one.tsv", rel)
            else:
                self.write(self.rsess, rel, rel)
        self.write(self.rsess, "output/search/report.pg_matrix.tsv")   # not on the list
        self.write(self.rsess, "logs/run.log")                          # not on the list

    def run_sync(self, *args, bash="bash", **env):
        e = {k: v for k, v in os.environ.items() if not k.startswith("HIVE_")}
        e.update(HIVE_EXEC=self.fake, FAKE_REMOTE=self.remote, FAKE_LOG=self.log, HOME=self.home,
                 TMPDIR=self.mk("tmp"))
        e.update(env)
        return subprocess.run([bash, SCRIPT, *args], capture_output=True, text=True, env=e,
                              timeout=120)

    def summary(self, r):
        return json.loads(r.stdout.strip().splitlines()[-1])

    def calls(self):
        if not os.path.exists(self.log):
            return []
        with open(self.log) as fh:
            return fh.read().splitlines()

    # ---- --get -------------------------------------------------------------------------
    def test_get_brings_every_item_to_its_place(self):
        self.full_session()
        r = self.run_sync("--get", "--session", self.SESSION, "--local", self.local)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        s = self.summary(r)
        self.assertEqual((s["mode"], s["session"], s["local"]), ("get", self.SESSION, self.local))
        self.assertEqual(s["done"], GET_ITEMS)
        self.assertEqual((s["missing"], s["failed"]), ([], []))
        for rel in GET_ITEMS:
            want = rel + "/one.tsv" if rel in ("output/tables", "output/figures") else rel
            with open(os.path.join(self.local, want)) as fh:
                self.assertEqual(fh.read(), rel, want)
            self.assertIn(f"got {self.SESSION}/{rel}\n", r.stdout)
        # nothing beyond the list, and nothing nested one level too deep
        self.assertFalse(os.path.exists(os.path.join(self.local, "output/search/report.pg_matrix.tsv")))
        self.assertFalse(os.path.exists(os.path.join(self.local, "logs")))
        self.assertFalse(os.path.exists(os.path.join(self.local, "output/tables/tables")))
        # one listing call, then one transfer per item
        self.assertEqual(len(self.calls()), 1 + len(GET_ITEMS))

    def test_absent_optional_items_are_missing_not_errors(self):
        self.full_session()
        for rel in ("output/AUDIT.md", "input/conditions.csv", "input/locate_decisions.json"):
            os.remove(os.path.join(self.rsess, rel))
        shutil.rmtree(os.path.join(self.rsess, "output", "figures"))
        r = self.run_sync("--get", "--session", self.SESSION, "--local", self.local)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        s = self.summary(r)
        gone = ["output/figures", "output/AUDIT.md", "input/conditions.csv",
                "input/locate_decisions.json"]
        self.assertEqual(s["missing"], gone)
        self.assertEqual(s["done"], [i for i in GET_ITEMS if i not in gone])
        self.assertEqual(s["failed"], [])
        for rel in gone:
            self.assertIn(f"missing {self.SESSION}/{rel}\n", r.stdout)
            self.assertFalse(os.path.exists(os.path.join(self.local, rel)))
            self.assertFalse(any(c.startswith("--get") and c.split()[1].endswith(rel)
                                 for c in self.calls()), rel)   # never asked for
        # the folders --get promises exist even when nothing landed in them
        self.assertTrue(os.path.isdir(os.path.join(self.local, "output", "search")))
        self.assertTrue(os.path.isdir(os.path.join(self.local, "input")))

    def test_an_absolute_session_path_and_crlf_from_hive(self):
        self.full_session()
        r = self.run_sync("--get", "--session", self.rsess + "/", "--local", self.local + "/",
                          FAKE_CRLF="1")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        s = self.summary(r)
        self.assertEqual((s["session"], s["local"]), (self.rsess, self.local))
        self.assertEqual(s["done"], GET_ITEMS)

    def test_the_session_records_are_in_the_get_list(self):
        """Moved from test_skill_md_session_records.py (it pinned SKILL.md's loop line): the
        records the local report reads come back with it."""
        r = self.run_sync("--get", "--dry-run", "--session", self.SESSION, "--local", self.local)
        self.assertEqual(r.returncode, 0, r.stderr)
        would = self.summary(r)["would"]
        for f in ("sample_type.json", "experiment_type.json", "locate_decisions.json"):
            self.assertIn("input/" + f, would)
        self.assertEqual(would, GET_ITEMS)

    def test_a_failed_transfer_is_failed_and_exits_1(self):
        self.full_session()
        r = self.run_sync("--get", "--session", self.SESSION, "--local", self.local,
                          FAKE_FAIL="output/figures")
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        s = self.summary(r)
        self.assertEqual(s["failed"], ["output/figures"])
        self.assertEqual(s["done"], [i for i in GET_ITEMS if i != "output/figures"])
        self.assertIn(f"failed {self.SESSION}/output/figures (hive_exec.sh --get exit 255)", r.stdout)

    def test_hive_unreachable_is_a_failure_not_everything_missing(self):
        self.full_session()
        r = self.run_sync("--get", "--session", self.SESSION, "--local", self.local, FAKE_DOWN="1")
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        s = self.summary(r)
        self.assertEqual((s["done"], s["missing"], s["failed"]), ([], [], ["."]))
        self.assertIn("exit 255", s["error"])
        self.assertEqual(len(self.calls()), 1)   # the listing call only

    def test_a_session_not_on_hive_fails(self):
        for mode in ("--get", "--put"):
            with self.subTest(mode=mode):
                self.write(self.local, "output/Analysis_Report.html")
                r = self.run_sync(mode, "--session", "~/core/PROT_0807/typo", "--local", self.local)
                self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
                s = self.summary(r)
                self.assertEqual((s["done"], s["failed"]), ([], ["."]))
                self.assertIn("no such session directory on HIVE", s["error"])
                # --put through scp would mkdir -p a stray tree under a typo: never got that far
                self.assertFalse(os.path.exists(os.path.join(self.remote, "core/PROT_0807/typo")))

    # ---- --put -------------------------------------------------------------------------
    def test_put_sends_only_existing_files_and_the_podcast(self):
        self.full_session()
        here = ["output/AI_Analysis_Report.md", "output/Analysis_Report.html", "output/methods.md",
                "output/methods.docx"]
        for rel in here:
            self.write(self.local, rel, "new " + rel)
        self.write(self.local, "output/podcast/podcast.m4a", "audio")
        self.write(self.local, "output/notes.txt")              # not a deliverable
        self.write(self.local, "output/tables/local_only.tsv")  # tables go HIVE -> here only
        r = self.run_sync("--put", "--session", self.SESSION, "--local", self.local)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        s = self.summary(r)
        self.assertEqual(s["done"], [i for i in PUT_ITEMS if i in here or i == "output/podcast"])
        self.assertEqual(s["missing"], ["output/Analysis_Report.md", "output/Analysis_Report.pdf",
                                        "output/Analysis_Report_with_audio.html"])
        self.assertEqual(s["failed"], [])
        for rel in here:
            with open(os.path.join(self.rsess, rel)) as fh:
                self.assertEqual(fh.read(), "new " + rel)
            self.assertIn(f"put {self.local}/{rel}\n", r.stdout)
        self.assertTrue(os.path.isfile(os.path.join(self.rsess, "output/podcast/podcast.m4a")))
        self.assertFalse(os.path.exists(os.path.join(self.rsess, "output/podcast/podcast")))
        self.assertIn(f"skipped {self.local}/output/Analysis_Report.pdf (not here)\n", r.stdout)
        self.assertFalse(os.path.exists(os.path.join(self.rsess, "output/notes.txt")))
        self.assertFalse(os.path.exists(os.path.join(self.rsess, "output/tables/local_only.tsv")))

    def test_put_without_an_output_folder_is_a_usage_error(self):
        os.makedirs(self.local)
        r = self.run_sync("--put", "--session", self.SESSION, "--local", self.local)
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertEqual(self.calls(), [])

    # ---- --dry-run, usage --------------------------------------------------------------
    def test_dry_run_transfers_nothing(self):
        self.full_session()
        r = self.run_sync("--get", "--dry-run", "--session", self.SESSION, "--local", self.local)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn(f"would get {self.SESSION}/input/conditions.csv -> {self.local}/input/", r.stdout)
        self.assertTrue(self.summary(r)["dry_run"])
        self.assertFalse(os.path.exists(self.local))
        self.write(self.local, "output/Analysis_Report.html")
        before = sorted(os.walk(self.remote))
        r = self.run_sync("--put", "--dry-run", "--session", self.SESSION, "--local", self.local)
        self.assertEqual(r.returncode, 0, r.stderr)
        s = self.summary(r)
        self.assertEqual((s["would"], s["done"]), (["output/Analysis_Report.html"], []))
        self.assertIn(f"would put {self.local}/output/Analysis_Report.html -> {self.SESSION}/output/",
                      r.stdout)
        self.assertEqual(sorted(os.walk(self.remote)), before)
        self.assertEqual(self.calls(), [])   # not even the listing call

    def test_bad_session_values_exit_2(self):
        bad = ["~/x;rm -rf ~", "/a/$(id)", "/a/`id`", '/a/"b', "/a/'b", "/a\\b", "/a|b", "/a&b",
               "/a\nb", "/a b", "/a/*", "relative/session", "~", "~/", "/", "~other/x", "/a/~/b", ""]
        for v in bad:
            with self.subTest(session=v):
                r = self.run_sync("--get", "--session", v, "--local", self.local)
                self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
                self.assertIn("usage:", r.stderr)
        self.assertEqual(self.calls(), [])
        self.assertFalse(os.path.exists(self.local))

    def test_usage_errors_exit_2(self):
        s, loc = ("--session", self.SESSION), ("--local", self.local)
        for args in ((*s, *loc), ("--get", *loc), ("--get", *s), ("--get", "--put", *s, *loc),
                     ("--get", *s, *loc, "--bogus"), ("--get", *s, "--local")):
            with self.subTest(args=args):
                r = self.run_sync(*args)
                self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertEqual(self.calls(), [])

    def test_help(self):
        for flag in ("-h", "--help"):
            r = self.run_sync(flag)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("usage: core_sync.sh --get|--put --session", r.stdout)
            self.assertIn("--dry-run", r.stdout)

    def test_local_tilde_and_json_escaping(self):
        self.full_session()
        r = self.run_sync("--get", "--session", self.SESSION, "--local", '~/a "q" \\ b')
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        want = os.path.join(self.home, 'a "q" \\ b')
        self.assertEqual(self.summary(r)["local"], want)
        self.assertTrue(os.path.isfile(os.path.join(want, "input", "sample_type.json")))

    def test_round_trip_under_every_bash_here(self):
        """macOS ships bash 3.2 as /bin/bash; Git Bash and CI are bash 4/5."""
        for n, b in enumerate(bashes()):
            with self.subTest(bash=b):
                self.layout(f"bash{n}")
                self.full_session()
                r = self.run_sync("--get", "--session", self.SESSION, "--local", self.local, bash=b)
                self.assertEqual((r.returncode, self.summary(r)["done"]), (0, GET_ITEMS), r.stderr)
                self.write(self.local, "output/Analysis_Report.html", "final")
                r = self.run_sync("--put", "--session", self.SESSION, "--local", self.local, bash=b)
                self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
                self.assertIn("output/Analysis_Report.html", self.summary(r)["done"])
                with open(os.path.join(self.rsess, "output/Analysis_Report.html")) as fh:
                    self.assertEqual(fh.read(), "final")

    # ---- the file itself ---------------------------------------------------------------
    def test_the_script_is_executable_and_parses(self):
        with open(SCRIPT) as fh:
            self.assertTrue(fh.readline().startswith("#!"))
        mode = os.stat(SCRIPT).st_mode
        self.assertTrue(mode & stat.S_IXUSR and mode & stat.S_IXGRP and mode & stat.S_IXOTH, oct(mode))
        p = subprocess.run(["bash", "-n", SCRIPT], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)


if __name__ == "__main__":
    unittest.main()
