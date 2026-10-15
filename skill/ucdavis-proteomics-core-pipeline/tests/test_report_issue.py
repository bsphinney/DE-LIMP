#!/usr/bin/env python3
"""
report_issue.sh is how a session's skill defects reach the people who fix them. The ten
defects of 2026-09-23 (a Windows Core member in hive_remote mode) were only fixed because the
agent happened to leave a notes file in the user's folder; this script makes that routine.

What must hold, because a failure here is silent -- nobody is waiting for the report:
  * three routes -- write to the Core's shared folder directly (on HIVE), through
    hive_exec.sh over SSH (a laptop in hive_remote mode), or into a local folder with a note
    to send it on -- and a failed shared write falls back to local instead of losing it;
  * one file per ENTRY, <date>_<user>[_<session>]+<unique>.md, written under a dot-name and
    renamed into place. Appending to one file per day lost entries to a check-then-create race
    (68 of 120 with 3 writers, release review 2.8.0), and a cross-node `>>` on /quobyte is
    unlocked -- flock does not lock across HIVE nodes there;
  * an entry without optional fields is still written (a `[ -n "$FIX" ] && printf` as the
    last line of the block once made every such entry look like a failed write);
  * text that looks like a key, token or password is refused, and nothing is written;
  * the session name cannot steer the file out of the folder.

Nothing here contacts HIVE: HIVE_EXEC is replaced by a stub that runs the remote command
against a temp directory standing in for /quobyte.
"""
import os
import re
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "scripts", "report_issue.sh")
FAKE_ROOT = "/nonexistent_hive_root_for_tests"     # never exists locally -> not "direct"


class Harness(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.shared = os.path.join(self.d, "shared")
        self.local = os.path.join(self.d, "local")
        self.remote = os.path.join(self.d, "remote")          # stands in for FAKE_ROOT on "HIVE"
        self.key = os.path.join(self.d, "id_test")
        open(self.key, "w").close()
        # hive_exec.sh stand-in: run the command locally with FAKE_ROOT mapped to self.remote,
        # stdin passed through exactly as ssh would.
        self.fake_exec = os.path.join(self.d, "fake_hive_exec.sh")
        with open(self.fake_exec, "w") as fh:
            fh.write('#!/usr/bin/env bash\ncmd="${1//%s/%s}"\nexec bash -c "$cmd"\n'
                     % (FAKE_ROOT, self.remote))
        os.chmod(self.fake_exec, 0o755)

    def tearDown(self):
        self._tmp.cleanup()

    def run_it(self, *args, route="direct", extra_env=None, stdin=None):
        env = {k: v for k, v in os.environ.items() if not k.startswith("HIVE_")}
        env.update(HOME=self.d, HIVE_ENV_FILE=os.path.join(self.d, "no-hive.env"),
                   SKILL_ISSUES_LOCAL_DIR=self.local, TMPDIR=self.d)
        if route == "direct":
            env["SKILL_ISSUES_DIR"] = self.shared
        else:
            env["SKILL_ISSUES_DIR"] = FAKE_ROOT + "/skill_issues"
        if route == "ssh":
            env.update(HIVE_USER="gabrig", HIVE_KEY=self.key, HIVE_EXEC=self.fake_exec)
        env.update(extra_env or {})
        return subprocess.run(["bash", SCRIPT, *args], capture_output=True, text=True, env=env,
                              input=stdin)

    def files(self, d):
        return sorted(os.listdir(d)) if os.path.isdir(d) else []

    def read_only_file(self, d):
        names = [n for n in self.files(d) if n.endswith(".md")]
        self.assertEqual(len(names), 1, names)
        with open(os.path.join(d, names[0])) as fh:
            return names[0], fh.read()


class ReportIssue(Harness):
    # ------------------------------------------------------------------ routes
    def test_each_entry_is_its_own_complete_file(self):
        os.makedirs(self.shared)
        for title in ("first", "second"):
            p = self.run_it("--title", title, "--what", "w", "--session", "chkLUppm")
            self.assertEqual(p.returncode, 0, p.stderr)
        names = self.files(self.shared)
        self.assertEqual(len(names), 2, names)
        titles = set()
        for n in names:
            # what the readers match on is the part before `+`
            self.assertRegex(n, r"^\d{4}-\d{2}-\d{2}_[^+]+_chkLUppm\+[A-Za-z0-9._-]+\.md$")
            with open(os.path.join(self.shared, n)) as fh:
                text = fh.read()
            self.assertEqual(text.count("\n## "), 1, text)
            self.assertIn("**Skill version:**", text)
            titles |= set(re.findall(r"^## (\w+)", text, re.M))
        self.assertEqual(titles, {"first", "second"})

    def test_concurrent_writers_lose_nothing(self):
        """Many writers at once, as sessions on several nodes are. The old one-file-per-day
        design raced at file creation -- every writer that found it empty truncated it: the
        release review lost 68 of 120, and 120 simultaneous writers here lost 14 of 120
        (2026-09-28). Every entry must land, whole, with nothing temporary left behind."""
        os.makedirs(self.shared)
        n = 60
        env = {k: v for k, v in os.environ.items() if not k.startswith("HIVE_")}
        env.update(HOME=self.d, HIVE_ENV_FILE=os.path.join(self.d, "no-hive.env"),
                   SKILL_ISSUES_LOCAL_DIR=self.local, TMPDIR=self.d, SKILL_ISSUES_DIR=self.shared)
        procs = [subprocess.Popen(["bash", SCRIPT, "--title", f"w{i}", "--what", "race",
                                   "--session", "race"], env=env, stdout=subprocess.DEVNULL,
                                  stderr=subprocess.PIPE, text=True) for i in range(n)]
        done = [(p, p.communicate()[1]) for p in procs]      # communicate() closes the pipes
        errors = [err for p, err in done if p.returncode != 0]
        self.assertEqual(errors, [])
        names = self.files(self.shared)
        self.assertEqual(len(names), n, "entries lost")
        self.assertFalse([x for x in names if x.startswith(".") or x.endswith(".part")])
        seen = set()
        for name in names:
            with open(os.path.join(self.shared, name)) as fh:
                seen |= set(re.findall(r"^## (w\d+)", fh.read(), re.M))
        self.assertEqual(len(seen), n)

    def test_entry_without_optional_fields_is_written(self):
        os.makedirs(self.shared)
        p = self.run_it("--title", "bare", "--what", "only the required fields")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertNotIn("cannot write", p.stderr)
        _, text = self.read_only_file(self.shared)
        self.assertIn("## bare", text)
        self.assertNotIn("Proposed fix", text)

    def test_ssh_route_writes_on_hive(self):
        os.makedirs(os.path.join(self.remote, "skill_issues"))
        p = self.run_it("--title", "over ssh", "--what", "w", "--fix", "f", route="ssh")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("recorded: gabrig@hive:", p.stdout)
        name, text = self.read_only_file(os.path.join(self.remote, "skill_issues"))
        self.assertTrue(name.startswith(tuple("0123456789")) and "_gabrig" in name, name)
        self.assertIn("hive_remote", text)
        self.assertEqual(self.files(self.local), [])

    def test_ssh_failure_falls_back_to_local(self):
        # no skill_issues folder on "HIVE" (a non-Core account): kept locally, not lost
        p = self.run_it("--title", "kept", "--what", "w", route="ssh")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("recorded locally", p.stdout)
        _, text = self.read_only_file(self.local)
        self.assertIn("## kept", text)

    def test_no_hive_at_all_goes_local_with_instructions(self):
        p = self.run_it("--title", "laptop", "--what", "w", route="local")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("send that file", p.stdout)
        self.read_only_file(self.local)

    def test_where_names_the_route(self):
        os.makedirs(self.shared)
        self.assertTrue(self.run_it("--where").stdout.startswith("direct:"))
        self.assertTrue(self.run_it("--where", route="ssh").stdout.startswith("ssh:"))
        self.assertTrue(self.run_it("--where", route="local").stdout.startswith("local:"))

    # ------------------------------------------------------------------ guards
    def test_secrets_are_refused_and_nothing_written(self):
        os.makedirs(self.shared)
        for bad in ("-----BEGIN OPENSSH PRIVATE KEY-----", "ghp_" + "a" * 30,
                    "Authorization: Token 0123456789abcdef", "password=hunter2"):
            p = self.run_it("--title", "t", "--what", bad)
            self.assertEqual(p.returncode, 2, bad)
            self.assertIn("key, token, password or webhook", p.stderr)
        self.assertEqual(self.files(self.shared), [])

    def test_session_name_cannot_escape_the_folder(self):
        os.makedirs(self.shared)
        p = self.run_it("--title", "t", "--what", "w", "--session", "../../etc/x y")
        self.assertEqual(p.returncode, 0, p.stderr)
        name, _ = self.read_only_file(self.shared)
        self.assertNotIn("/", name)
        self.assertNotIn(" ", name)

    def test_required_fields_and_enums(self):
        self.assertEqual(self.run_it("--title", "t").returncode, 2)
        self.assertEqual(self.run_it("--title", "t", "--what", "w", "--kind", "oops").returncode, 2)
        self.assertEqual(self.run_it("--title", "t", "--what", "w",
                                     "--severity", "urgent").returncode, 2)


class TextFromAFile(Harness):
    """2.12 security review: the exact error text went on the command line (--what "..."), where
    bash runs anything in backticks or $(...) -- on this computer, and again on HIVE through
    hive_exec.sh's `bash -l -c` (a backtick ran a command). The text now comes from a file or
    standard input (--what-file F|-, and the same for the other text fields), read as data."""

    def hostile(self):
        """An error text that runs commands if a shell ever evaluates it."""
        self.ran = [os.path.join(self.d, n) for n in ("RAN_BACKTICK", "RAN_DOLLAR")]
        return (f'ls failed: `touch {self.ran[0]}` and $(touch {self.ran[1]})\n'
                f'second line: "double" \'single\' $HOME ${{PATH}} \\ end\n')

    def entry_field(self, text, label):
        return text.split(f"- **{label}:** ", 1)[1].split("\n- **", 1)[0].rstrip("\n")

    def test_the_exact_text_from_stdin_is_written_and_never_run(self):
        os.makedirs(self.shared)
        what = self.hostile()
        p = self.run_it("--title", "t", "--what-file", "-", "--impact", "i", stdin=what)
        self.assertEqual(p.returncode, 0, p.stderr)
        for f in self.ran:
            self.assertFalse(os.path.exists(f), "the error text was run as a command")
        _name, text = self.read_only_file(self.shared)
        lines = what.rstrip("\n").split("\n")
        self.assertEqual(self.entry_field(text, "What happened"),
                         lines[0] + "\n  " + lines[1])     # verbatim, under its bullet

    def test_over_ssh_too(self):
        os.makedirs(os.path.join(self.remote, "skill_issues"))
        what = self.hostile()
        p = self.run_it("--title", "t", "--what-file", "-", stdin=what, route="ssh")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("recorded: gabrig@hive", p.stdout)
        for f in self.ran:
            self.assertFalse(os.path.exists(f))
        _name, text = self.read_only_file(os.path.join(self.remote, "skill_issues"))
        self.assertIn(what.split("\n")[0], text)

    def test_every_text_field_can_come_from_a_file(self):
        os.makedirs(self.shared)
        paths = {}
        for field in ("title", "what", "impact", "workaround", "fix"):
            paths[field] = os.path.join(self.d, field + ".txt")
            with open(paths[field], "w") as fh:
                fh.write(f"the {field} `x` $(y)\r\nline two\n")
        p = self.run_it(*[a for f, path in paths.items() for a in (f"--{f}-file", path)])
        self.assertEqual(p.returncode, 0, p.stderr)
        _name, text = self.read_only_file(self.shared)
        self.assertIn("## the title `x` $(y) line two  [bug, medium]", text)    # a heading: one line
        for field, label in (("what", "What happened"), ("impact", "Impact"),
                             ("workaround", "Workaround used"), ("fix", "Proposed fix")):
            self.assertEqual(self.entry_field(text, label), f"the {field} `x` $(y)\n  line two")
        self.assertNotIn("\r", text)

    def test_secrets_in_a_file_are_refused_too(self):
        os.makedirs(self.shared)
        for bad in ("ghp_" + "a" * 30, "line one\nAuthorization: Token 0123456789abcdef\n"):
            p = self.run_it("--title", "t", "--what-file", "-", stdin=bad)
            self.assertEqual(p.returncode, 2, bad)
            self.assertIn("key, token, password or webhook", p.stderr)
        self.assertEqual(self.files(self.shared), [])

    def test_wrong_uses_are_refused(self):
        os.makedirs(self.shared)
        f = os.path.join(self.d, "w.txt")
        with open(f, "w") as fh:
            fh.write("w")
        for args, said in ((("--what-file", "-", "--fix-file", "-"), "only one --...-file can be -"),
                           (("--what", "w", "--what-file", f), "give --what or --what-file, not both"),
                           (("--what-file", os.path.join(self.d, "nope.txt")), "not a readable file"),
                           (("--what-file", self.d), "not a readable file")):
            with self.subTest(args=args):
                p = self.run_it("--title", "t", *args, stdin="x")
                self.assertEqual(p.returncode, 2, p.stdout)
                self.assertIn(said, p.stderr)
        p = self.run_it("--title", "t", "--what-file", "-", stdin="")
        self.assertEqual(p.returncode, 2)
        self.assertIn("are required", p.stderr)
        self.assertEqual(self.files(self.shared), [])

    def test_a_huge_text_is_cut_at_64_kib(self):
        os.makedirs(self.shared)
        p = self.run_it("--title", "t", "--what-file", "-", stdin="a" * 70000)
        self.assertEqual(p.returncode, 0, p.stderr)
        _name, text = self.read_only_file(self.shared)
        self.assertEqual(self.entry_field(text, "What happened"), "a" * 65536 + " [... cut at 64 KiB]")

    def test_help_and_skill_md_give_the_stdin_route(self):
        h = self.run_it("--help").stdout
        self.assertIn("--what-file -", h)
        self.assertIn("<<'EOF'", h)
        self.assertTrue(h.rstrip().endswith("# " + "=" * 77), h[-200:])
        skill = os.path.join(os.path.dirname(HERE), "SKILL.md")
        with open(skill, encoding="utf-8") as fh:
            doc = fh.read()
        sec = doc[doc.index("## Recording skill problems"):]
        sec = sec[:sec.index("\n## ", 5)]
        block = sec.split("```")[1]                      # the command shown
        self.assertIn("--what-file -", block)
        self.assertIn("<<'EOF'", block)
        self.assertNotRegex(block, r'--(title|what|impact|workaround|fix|step|session) "',
                            "a report_issue.sh text argument in double quotes")
        self.assertIn("never inside `--what \"...\"`", sec)


if __name__ == "__main__":
    unittest.main()
